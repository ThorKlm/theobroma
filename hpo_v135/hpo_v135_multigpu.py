"""v1.35 classifier HPO, multi-GPU (trial-parallel) with SQLite-backed shared study.

Parallelism model: several independent worker processes, each pinned to one GPU via
CUDA_VISIBLE_DEVICES, all pulling trials from ONE shared Optuna study (SQLite). Trials
are independent, so this is embarrassingly parallel; the TPE sampler coordinates via
the shared storage. No single fit is split across GPUs (unnecessary for this model).

Two roles (selected by --role):
  worker      : runs a slice of the trial budget, then exits. Launch one per GPU.
  coordinate  : waits until the study has >= --trials completed trials, then does the
                one-time finalization (repeated CV, holdout eval, baselines, model save).
                Run exactly one of these. It does NOT train trials itself.
A single-process convenience role 'all' reproduces the original single-GPU behaviour.

SQLite robustness: the storage is opened with a busy timeout and SQLAlchemy
pool_pre_ping, and trial writes are retried on 'database is locked'. With <=4 workers
and early-stopped fits (tens of seconds to minutes each), write contention is rare
because writes happen only at trial start/end, not during boosting.

Outputs (written once, by the coordinator): identical to the single-GPU script,
xgb_v135.ubj, pca_components.npy, pca_mean.npy, classes.json, hpo_results_v135.json,
hpo_per_class_eval_v135.csv, threshold_calibration_v135.csv, hpo_trial_log_v135.csv.
"""
import argparse, json, os, time
import numpy as np, pandas as pd
import optuna, xgboost as xgb
from scipy import sparse
from sklearn.decomposition import IncrementalPCA
from sklearn.metrics import f1_score, log_loss
from sklearn.model_selection import StratifiedKFold, train_test_split
from tqdm.auto import tqdm

# --- storage with SQLite busy-timeout + locked-write retry -------------------

def make_storage(path, busy_timeout_ms=60000):
    """Optuna RDBStorage on SQLite with a busy timeout and resilient engine kwargs.

    busy_timeout lets SQLite block (not error) for up to the timeout when another
    worker holds the write lock; heartbeat/grace let stalled trials be reclaimed.
    """
    return optuna.storages.RDBStorage(
        url=f"sqlite:///{path}",
        engine_kwargs={"connect_args": {"timeout": busy_timeout_ms / 1000.0},
                       "pool_pre_ping": True},
        heartbeat_interval=60, grace_period=180,
        failed_trial_callback=optuna.storages.RetryFailedTrialCallback(max_retry=3),
    )

def with_retry(fn, tries=8, base=0.5):
    """Retry a storage-touching call on transient 'database is locked'."""
    import sqlite3
    from sqlalchemy.exc import OperationalError
    for i in range(tries):
        try:
            return fn()
        except (OperationalError, sqlite3.OperationalError) as e:
            if "locked" not in str(e).lower() or i == tries - 1:
                raise
            time.sleep(base * (2 ** i))

# --- feature assembly -------------------------------------------------------

def build_features(vdir, row_idx, pca_dim=256, pca_sample=200_000, seed=0):
    morgan = sparse.load_npz(os.path.join(vdir, "morgan_fps.npz")).tocsr()[row_idx]
    maccs = sparse.load_npz(os.path.join(vdir, "maccs_fps.npz")).tocsr()[row_idx]
    emb = np.load(os.path.join(vdir, "chemberta_embeddings.npy"), mmap_mode="r")
    emb = np.asarray(emb[row_idx], dtype=np.float32)
    rng = np.random.default_rng(seed)
    fit_idx = rng.choice(len(emb), size=min(pca_sample, len(emb)), replace=False)
    pca = IncrementalPCA(n_components=pca_dim, batch_size=8192)
    pca.fit(emb[np.sort(fit_idx)])
    proj = pca.transform(emb).astype(np.float32)
    print(f"[feat] pca explained variance {pca.explained_variance_ratio_.sum():.3f}", flush=True)
    X = np.hstack([morgan.toarray().astype(np.float32), maccs.toarray().astype(np.float32), proj])
    return X, pca

def class_weights(y, mode):
    if mode == "none":
        return None
    counts = np.bincount(y)
    inv = counts.sum() / np.maximum(counts, 1)
    w = inv if mode == "full" else np.sqrt(inv)
    w = w / w.mean()
    return w[y].astype(np.float32)

# --- tqdm callback ----------------------------------------------------------

class TqdmCallback(xgb.callback.TrainingCallback):
    def __init__(self, total, desc, position=0):
        self.bar = tqdm(total=total, desc=desc, leave=False, dynamic_ncols=True,
                        unit="round", mininterval=0.5, position=position)
    def after_iteration(self, model, epoch, evals_log):
        msg = ""
        if evals_log:
            ds = next(iter(evals_log)); metric = next(iter(evals_log[ds]))
            msg = f"{ds}-{metric}={evals_log[ds][metric][-1]:.4f}"
        self.bar.set_postfix_str(msg); self.bar.update(1)
        return False
    def after_training(self, model):
        self.bar.close(); return model

# --- training ---------------------------------------------------------------

def _params(p, n_class, seed):
    p = dict(p); p.pop("weight_mode", None)
    p.update(objective="multi:softprob", num_class=n_class, eval_metric="mlogloss",
             tree_method="hist", device="cuda", seed=seed, max_bin=256)
    return p

def fit_on(dtr, dva, params, n_class, max_rounds, esr, seed, desc, position=0):
    p = _params(params, n_class, seed)
    bst = xgb.train(p, dtr, num_boost_round=max_rounds, evals=[(dva, "va")],
                    early_stopping_rounds=esr, verbose_eval=False,
                    callbacks=[TqdmCallback(max_rounds, desc, position)])
    proba = bst.predict(dva, iteration_range=(0, bst.best_iteration + 1))
    return bst, proba, bst.best_iteration + 1

def make_dmatrix(X, y, wmode, ref=None):
    return xgb.QuantileDMatrix(X, label=y, weight=class_weights(y, wmode), ref=ref)

def suggest(trial):
    return {
        "max_depth": trial.suggest_int("max_depth", 4, 14),
        "learning_rate": trial.suggest_float("learning_rate", 0.05, 0.5, log=True),
        "subsample": trial.suggest_float("subsample", 0.5, 1.0),
        "colsample_bytree": trial.suggest_float("colsample_bytree", 0.1, 1.0),
        "min_child_weight": trial.suggest_int("min_child_weight", 1, 30),
        "reg_lambda": trial.suggest_float("reg_lambda", 1e-2, 50.0, log=True),
        "reg_alpha": trial.suggest_float("reg_alpha", 1e-3, 20.0, log=True),
        "gamma": trial.suggest_float("gamma", 0.0, 5.0),
        "weight_mode": trial.suggest_categorical("weight_mode", ["none", "sqrt", "full"]),
    }

# --- evaluation -------------------------------------------------------------

def per_class_report(y_true, y_pred, classes):
    rows = []
    for k, name in enumerate(classes):
        tp = int(((y_pred == k) & (y_true == k)).sum())
        fp = int(((y_pred == k) & (y_true != k)).sum())
        fn = int(((y_pred != k) & (y_true == k)).sum())
        prec = tp / (tp + fp) if tp + fp else 0.0
        rec = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
        rows.append({"class": name, "support": int((y_true == k).sum()),
                     "precision": round(prec, 4), "recall": round(rec, 4), "f1": round(f1, 4)})
    return pd.DataFrame(rows)

def calibrate_thresholds(proba, y_true, grid=np.arange(0.05, 0.96, 0.05)):
    top = proba.max(axis=1); pred = proba.argmax(axis=1); rows = []
    for t in grid:
        m = top >= t
        rows.append({"threshold": round(float(t), 2), "coverage": round(float(m.mean()), 4),
                     "accuracy_on_accepted": round(float((pred[m] == y_true[m]).mean()), 4) if m.any() else None,
                     "macro_f1_on_accepted": round(float(f1_score(y_true[m], pred[m], average="macro",
                                                     zero_division=0)), 4) if m.any() else None})
    return pd.DataFrame(rows)

def cv_evaluate(X, y, params, n_class, classes, max_rounds, esr, seed, n_splits=5):
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    macro, micro, oof = [], [], np.zeros((len(y), n_class), dtype=np.float32)
    wmode = params.get("weight_mode", "none")
    for i, (tr, va) in enumerate(skf.split(X, y)):
        dtr = make_dmatrix(X[tr], y[tr], wmode)
        dva = make_dmatrix(X[va], y[va], "none", ref=dtr)
        _, proba, _ = fit_on(dtr, dva, params, n_class, max_rounds, esr, seed, f"cv {i+1}/{n_splits}")
        oof[va] = proba
        pred = proba.argmax(axis=1)
        macro.append(f1_score(y[va], pred, average="macro", zero_division=0))
        micro.append(f1_score(y[va], pred, average="micro", zero_division=0))
    pred = oof.argmax(axis=1)
    return {"macro_f1_mean": float(np.mean(macro)), "macro_f1_std": float(np.std(macro)),
            "micro_f1_mean": float(np.mean(micro)), "micro_f1_std": float(np.std(micro)),
            "folds": [round(float(m), 4) for m in macro]}, oof, pred

# --- data loading (shared) --------------------------------------------------

def load_split(args):
    lab = pd.read_csv(args.labels, sep="\t")
    if args.stratum == "native":
        lab = lab[lab.origin == "native"]
    elif args.stratum == "consistent":
        lab = lab[(lab.origin == "native") | (lab.pathway_consistent == 1)]
    counts = lab.np_class.value_counts()
    keep = counts[counts >= args.min_class].index
    lab = lab[lab.np_class.isin(keep)].reset_index(drop=True)
    classes = sorted(lab.np_class.unique())
    cmap = {c: i for i, c in enumerate(classes)}
    y = lab.np_class.map(cmap).to_numpy()
    print(f"[data] stratum={args.stratum} n={len(lab):,} classes={len(classes)} "
          f"(min support {counts[keep].min()}, max {counts[keep].max()})", flush=True)
    X, pca = build_features(args.vectors, lab.row_idx.to_numpy(), seed=args.seed)
    print(f"[data] X {X.shape} {X.nbytes / 1e9:.2f} GB", flush=True)
    Xs, Xh, ys, yh = train_test_split(X, y, test_size=0.15, stratify=y, random_state=args.seed)
    Xtr, Xva, ytr, yva = train_test_split(Xs, ys, test_size=0.18, stratify=ys, random_state=args.seed)
    return dict(lab=lab, classes=classes, X=X, pca=pca, Xs=Xs, Xh=Xh, ys=ys, yh=yh,
                Xtr=Xtr, Xva=Xva, ytr=ytr, yva=yva)

ENQUEUE = {"max_depth": 9, "learning_rate": 0.21, "subsample": 1.0, "colsample_bytree": 0.3,
           "min_child_weight": 2, "reg_lambda": 6.8, "reg_alpha": 6.7, "gamma": 0.0,
           "weight_mode": "sqrt"}

# --- worker role ------------------------------------------------------------

def run_worker(args, storage):
    d = load_split(args)
    classes, n_class = d["classes"], len(d["classes"])
    Xtr, Xva, ytr, yva = d["Xtr"], d["Xva"], d["ytr"], d["yva"]
    dtr_cache, dva_cache = {}, {}
    def get_split(wmode):
        if wmode not in dtr_cache:
            dtr = xgb.QuantileDMatrix(Xtr, label=ytr, weight=class_weights(ytr, wmode))
            dtr_cache[wmode] = dtr; dva_cache[wmode] = xgb.QuantileDMatrix(Xva, label=yva, ref=dtr)
        return dtr_cache[wmode], dva_cache[wmode]

    def objective(trial):
        p = suggest(trial); t0 = time.time()
        dtr, dva = get_split(p["weight_mode"])
        _, proba, rounds = fit_on(dtr, dva, p, n_class, args.max_rounds, args.early_stopping,
                                  args.seed, f"w{args.worker_id} t{trial.number} "
                                  f"d{p['max_depth']} lr{p['learning_rate']:.2g} {p['weight_mode']}")
        score = f1_score(yva, proba.argmax(axis=1), average="macro", zero_division=0)
        trial.set_user_attr("rounds", rounds)
        trial.set_user_attr("fit_min", round((time.time() - t0) / 60, 2))
        tqdm.write(f"[w{args.worker_id} trial {trial.number}] macroF1={score:.4f} "
                   f"rounds={rounds} fit={trial.user_attrs['fit_min']}min wmode={p['weight_mode']}")
        return score

    sampler = optuna.samplers.TPESampler(multivariate=True, group=True, seed=args.seed + args.worker_id,
                                         n_startup_trials=args.startup)
    study = optuna.create_study(direction="maximize", sampler=sampler, storage=storage,
                                study_name=args.study, load_if_exists=True)
    # each worker runs at most its slice; stop early if the global budget is met
    def budget_reached(_study, _trial):
        n_done = len([t for t in _study.trials if t.state.name in ("COMPLETE", "PRUNED")])
        if n_done >= args.trials:
            _study.stop()
    per_worker = max(1, args.trials // max(1, args.workers)) + 2  # small overshoot; global cap enforced
    with_retry(lambda: study.optimize(objective, n_trials=per_worker, callbacks=[budget_reached],
                                      gc_after_trial=True))
    print(f"[w{args.worker_id}] done", flush=True)

# --- coordinator role -------------------------------------------------------

def run_coordinate(args, storage):
    study = optuna.create_study(direction="maximize", storage=storage,
                                study_name=args.study, load_if_exists=True)
    print(f"[coord] waiting for {args.trials} completed trials ...", flush=True)
    while True:
        n = len([t for t in study.trials if t.state.name == "COMPLETE"])
        if n >= args.trials:
            break
        time.sleep(30)
    print(f"[coord] {args.trials} trials complete; finalizing", flush=True)
    d = load_split(args)  # coordinator rebuilds features (deterministic, same seed)
    classes, n_class = d["classes"], len(d["classes"])
    Xs, Xh, ys, yh, pca = d["Xs"], d["Xh"], d["ys"], d["yh"], d["pca"]
    study.trials_dataframe().to_csv(os.path.join(args.out, "hpo_trial_log_v135.csv"), index=False)
    best = dict(study.best_params)
    print(f"[hpo ] best validation macro-F1 {study.best_value:.4f}", flush=True)
    cv, oof, oof_pred = cv_evaluate(Xs, ys, best, n_class, classes,
                                    args.max_rounds, args.early_stopping, args.seed)
    print(f"[cv  ] macro-F1 {cv['macro_f1_mean']:.4f} +/- {cv['macro_f1_std']:.4f}", flush=True)
    dtr_h = make_dmatrix(Xs, ys, best.get("weight_mode", "none"))
    dva_h = make_dmatrix(Xh, yh, "none", ref=dtr_h)
    bst, hproba, rounds = fit_on(dtr_h, dva_h, best, n_class, args.max_rounds,
                                 args.early_stopping, args.seed, "holdout fit")
    hpred = hproba.argmax(axis=1)
    per_class_report(yh, hpred, classes).to_csv(
        os.path.join(args.out, "hpo_per_class_eval_v135.csv"), index=False)
    calibrate_thresholds(hproba, yh).to_csv(
        os.path.join(args.out, "threshold_calibration_v135.csv"), index=False)
    maj = np.bincount(ys, minlength=n_class).argmax()
    baselines = {"majority_class_macro_f1": float(f1_score(yh, np.full_like(yh, maj),
                                                  average="macro", zero_division=0))}
    for wm in ("none", "sqrt", "full"):
        p = dict(best, weight_mode=wm)
        d1 = make_dmatrix(Xs, ys, wm); d2 = make_dmatrix(Xh, yh, "none", ref=d1)
        _, pr, _ = fit_on(d1, d2, p, n_class, args.max_rounds, args.early_stopping, args.seed, f"baseline {wm}")
        baselines[f"weight_mode_{wm}_macro_f1"] = float(f1_score(yh, pr.argmax(axis=1),
                                                        average="macro", zero_division=0))
    results = {
        "config": {"stratum": args.stratum, "min_class": args.min_class, "n_classes": n_class,
                   "n_labelled": int(len(d["lab"])), "n_search": int(len(d["Xtr"])),
                   "n_holdout": int(len(Xh)), "trials": args.trials,
                   "feature_dim": int(d["X"].shape[1]), "seed": args.seed, "workers": args.workers},
        "best": {"params": best, "rounds": rounds, "search_macro_f1": float(study.best_value)},
        "cv": cv,
        "holdout": {"macro_f1": float(f1_score(yh, hpred, average="macro", zero_division=0)),
                    "micro_f1": float(f1_score(yh, hpred, average="micro", zero_division=0)),
                    "logloss": float(log_loss(yh, hproba, labels=list(range(n_class))))},
        "baselines": baselines,
    }
    with open(os.path.join(args.out, "hpo_results_v135.json"), "w") as fh:
        json.dump(results, fh, indent=1)
    bst.save_model(os.path.join(args.out, "xgb_v135.ubj"))
    np.save(os.path.join(args.out, "pca_components.npy"), pca.components_)
    np.save(os.path.join(args.out, "pca_mean.npy"), pca.mean_)
    with open(os.path.join(args.out, "classes.json"), "w") as fh:
        json.dump(classes, fh)
    print(json.dumps(results["holdout"], indent=1), flush=True)

# --- main -------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vectors", default="data/vectors")
    ap.add_argument("--labels", default="data/vectors/labels_v135.tsv")
    ap.add_argument("--out", default="out_v135")
    ap.add_argument("--stratum", choices=["native", "all", "consistent"], default="native")
    ap.add_argument("--min-class", type=int, default=50)
    ap.add_argument("--trials", type=int, default=75, help="GLOBAL completed-trial budget (startup+TPE)")
    ap.add_argument("--startup", type=int, default=50, help="random exploration trials (pooled)")
    ap.add_argument("--workers", type=int, default=4, help="number of GPU workers")
    ap.add_argument("--worker-id", type=int, default=0)
    ap.add_argument("--role", choices=["worker", "coordinate", "all"], default="all")
    ap.add_argument("--max-rounds", type=int, default=500)
    ap.add_argument("--early-stopping", type=int, default=25)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--storage-path", default="optuna_v135.db")
    ap.add_argument("--study", default="theobroma_v135")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    storage = make_storage(args.storage_path)
    if args.role == "worker":
        run_worker(args, storage)
    elif args.role == "coordinate":
        run_coordinate(args, storage)
    else:
        # single-process: optimize then finalize (original behaviour)
        args.workers = 1
        run_worker(args, storage)
        run_coordinate(args, storage)

if __name__ == "__main__":
    main()
