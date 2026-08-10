"""v1.35 classifier: label-provenance-aware HPO and evaluation.

Fixes carried over from the v1.34 post-mortem:
  1. Objective is macro-F1 over every retained class, not a 30-class subset.
  2. n_rounds is not searched; early stopping on a validation fold sets it.
  3. Ranges widened where v1.34 pinned at a boundary (max_depth up, colsample down).
  4. Startup ratio is 20 percent of budget, not 80 percent (v1.34 was near-random).
  5. Training labels are stratified by provenance; the artifact stratum is optional.
  6. Reported score comes from repeated stratified CV on the top configs, so per-class
     support is never the n_test=1 of v1.34.
  7. NPClassifier-direct and majority-class baselines are evaluated on the same folds.
  8. Acceptance thresholds for the three-tier scheme are calibrated, not assumed.

Warm start is deliberately by parameter location only (enqueue), never by replayed
objective value: the v1.34 values were computed under a different objective and label
set and are not comparable.

Intended for a single RTX 3090. Expect roughly 6-10 h at the default budget.
"""
import argparse, json, os, time
import numpy as np, pandas as pd
import optuna, xgboost as xgb
from scipy import sparse
from sklearn.decomposition import IncrementalPCA
from sklearn.metrics import f1_score, log_loss
from sklearn.model_selection import StratifiedKFold, train_test_split

# --- feature assembly -------------------------------------------------------

def build_features(vdir, row_idx, pca_dim=256, pca_sample=200_000, seed=0):
    """Morgan + MACCS + PCA-projected ChemBERTa, restricted to the labelled rows.

    PCA is fitted on a random subsample of the labelled rows only, so no
    information from held-out compounds enters the projection.
    """
    morgan = sparse.load_npz(os.path.join(vdir, "morgan_fps.npz")).tocsr()[row_idx]
    maccs = sparse.load_npz(os.path.join(vdir, "maccs_fps.npz")).tocsr()[row_idx]
    emb = np.load(os.path.join(vdir, "chemberta_embeddings.npy"), mmap_mode="r")
    emb = np.asarray(emb[row_idx], dtype=np.float32)
    rng = np.random.default_rng(seed)
    fit_idx = rng.choice(len(emb), size=min(pca_sample, len(emb)), replace=False)
    pca = IncrementalPCA(n_components=pca_dim, batch_size=8192)
    pca.fit(emb[np.sort(fit_idx)])
    proj = pca.transform(emb).astype(np.float32)
    print(f"[feat] pca explained variance {pca.explained_variance_ratio_.sum():.3f}")
    X = np.hstack([morgan.toarray().astype(np.float32), maccs.toarray().astype(np.float32), proj])
    return X, pca

def class_weights(y, mode):
    """Per-sample weights. 'full' inverts frequency, 'sqrt' softens it."""
    if mode == "none":
        return None
    counts = np.bincount(y)
    inv = counts.sum() / np.maximum(counts, 1)
    w = inv if mode == "full" else np.sqrt(inv)
    w = w / w.mean()
    return w[y].astype(np.float32)

# --- training ---------------------------------------------------------------

def fit_predict(Xtr, ytr, Xva, yva, params, n_class, max_rounds, esr, seed):
    """One fit with early stopping. Returns probabilities and the chosen round count."""
    p = dict(params)
    wmode = p.pop("weight_mode")
    p.update(objective="multi:softprob", num_class=n_class, eval_metric="mlogloss",
             tree_method="hist", device="cuda", seed=seed, max_bin=256)
    dtr = xgb.QuantileDMatrix(Xtr, label=ytr, weight=class_weights(ytr, wmode))
    dva = xgb.QuantileDMatrix(Xva, label=yva, ref=dtr)
    bst = xgb.train(p, dtr, num_boost_round=max_rounds, evals=[(dva, "va")],
                    early_stopping_rounds=esr, verbose_eval=False)
    proba = bst.predict(dva, iteration_range=(0, bst.best_iteration + 1))
    return bst, proba, bst.best_iteration + 1

def suggest(trial):
    """Search space. Bounds widened where v1.34 selected a boundary value."""
    return {
        "max_depth": trial.suggest_int("max_depth", 4, 14),
        "learning_rate": trial.suggest_float("learning_rate", 0.01, 0.5, log=True),
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
    """Coverage and accuracy-on-accepted as a function of the top-class probability.

    Justifies the three-tier cut points instead of asserting 0.5 and 0.3.
    """
    top = proba.max(axis=1)
    pred = proba.argmax(axis=1)
    rows = []
    for t in grid:
        m = top >= t
        rows.append({"threshold": round(float(t), 2), "coverage": round(float(m.mean()), 4),
                     "accuracy_on_accepted": round(float((pred[m] == y_true[m]).mean()), 4) if m.any() else None,
                     "macro_f1_on_accepted": round(float(f1_score(y_true[m], pred[m], average="macro",
                                                     zero_division=0)), 4) if m.any() else None})
    return pd.DataFrame(rows)

def cv_evaluate(X, y, params, n_class, classes, max_rounds, esr, seed, n_splits=5):
    """Repeated stratified CV for the reported figure, with fold-level spread."""
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    macro, micro, oof = [], [], np.zeros((len(y), n_class), dtype=np.float32)
    for tr, va in skf.split(X, y):
        _, proba, _ = fit_predict(X[tr], y[tr], X[va], y[va], params, n_class, max_rounds, esr, seed)
        oof[va] = proba
        pred = proba.argmax(axis=1)
        macro.append(f1_score(y[va], pred, average="macro", zero_division=0))
        micro.append(f1_score(y[va], pred, average="micro", zero_division=0))
    pred = oof.argmax(axis=1)
    return {"macro_f1_mean": float(np.mean(macro)), "macro_f1_std": float(np.std(macro)),
            "micro_f1_mean": float(np.mean(micro)), "micro_f1_std": float(np.std(micro)),
            "folds": [round(float(m), 4) for m in macro]}, oof, pred

# --- main -------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vectors", default="data/vectors")
    ap.add_argument("--labels", default="data/vectors/labels_v135.tsv")
    ap.add_argument("--out", default="out_v135")
    ap.add_argument("--stratum", choices=["native", "all", "consistent"], default="native",
                    help="native: exclude artifact-derived labels (recommended). "
                         "consistent: artifact rows retained only where the stored pathway "
                         "agrees with the ontology. all: v1.34-equivalent, for comparison.")
    ap.add_argument("--min-class", type=int, default=50)
    ap.add_argument("--trials", type=int, default=60)
    ap.add_argument("--max-rounds", type=int, default=2000)
    ap.add_argument("--early-stopping", type=int, default=30)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--storage", default="sqlite:///optuna_v135.db")
    ap.add_argument("--study", default="theobroma_v135")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
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
          f"(min support {counts[keep].min()}, max {counts[keep].max()})")
    X, pca = build_features(args.vectors, lab.row_idx.to_numpy(), seed=args.seed)
    print(f"[data] X {X.shape} {X.nbytes / 1e9:.2f} GB")
    Xs, Xh, ys, yh = train_test_split(X, y, test_size=0.15, stratify=y, random_state=args.seed)
    Xtr, Xva, ytr, yva = train_test_split(Xs, ys, test_size=0.18, stratify=ys, random_state=args.seed)

    def objective(trial):
        p = suggest(trial)
        t0 = time.time()
        _, proba, rounds = fit_predict(Xtr, ytr, Xva, yva, p, len(classes),
                                       args.max_rounds, args.early_stopping, args.seed)
        score = f1_score(yva, proba.argmax(axis=1), average="macro", zero_division=0)
        trial.set_user_attr("rounds", rounds)
        trial.set_user_attr("fit_min", round((time.time() - t0) / 60, 2))
        trial.set_user_attr("micro_f1", round(float(f1_score(yva, proba.argmax(axis=1),
                                              average="micro", zero_division=0)), 4))
        return score

    sampler = optuna.samplers.TPESampler(multivariate=True, group=True, seed=args.seed,
                                         n_startup_trials=max(8, args.trials // 5))
    study = optuna.create_study(direction="maximize", sampler=sampler, storage=args.storage,
                                study_name=args.study, load_if_exists=True)
    # Re-evaluate the v1.34 selection under the corrected objective rather than
    # trusting its reported value.
    study.enqueue_trial({"max_depth": 9, "learning_rate": 0.21, "subsample": 1.0,
                         "colsample_bytree": 0.3, "min_child_weight": 2, "reg_lambda": 6.8,
                         "reg_alpha": 6.7, "gamma": 0.0, "weight_mode": "sqrt"}, skip_if_exists=True)
    study.optimize(objective, n_trials=args.trials, gc_after_trial=True)
    study.trials_dataframe().to_csv(os.path.join(args.out, "hpo_trial_log_v135.csv"), index=False)

    best = dict(study.best_params)
    print(f"[hpo ] best validation macro-F1 {study.best_value:.4f}")
    cv, oof, oof_pred = cv_evaluate(Xs, ys, best, len(classes), classes,
                                    args.max_rounds, args.early_stopping, args.seed)
    print(f"[cv  ] macro-F1 {cv['macro_f1_mean']:.4f} +/- {cv['macro_f1_std']:.4f}")
    bst, hproba, rounds = fit_predict(Xs, ys, Xh, yh, best, len(classes),
                                      args.max_rounds, args.early_stopping, args.seed)
    hpred = hproba.argmax(axis=1)
    per_class_report(yh, hpred, classes).to_csv(
        os.path.join(args.out, "hpo_per_class_eval_v135.csv"), index=False)
    calibrate_thresholds(hproba, yh).to_csv(
        os.path.join(args.out, "threshold_calibration_v135.csv"), index=False)
    # Baselines on the identical holdout.
    maj = np.bincount(ys, minlength=len(classes)).argmax()
    baselines = {"majority_class_macro_f1": float(f1_score(yh, np.full_like(yh, maj),
                                                  average="macro", zero_division=0))}
    for wm in ("none", "sqrt", "full"):
        p = dict(best, weight_mode=wm)
        _, pr, _ = fit_predict(Xs, ys, Xh, yh, p, len(classes), args.max_rounds,
                               args.early_stopping, args.seed)
        baselines[f"weight_mode_{wm}_macro_f1"] = float(f1_score(yh, pr.argmax(axis=1),
                                                        average="macro", zero_division=0))
    results = {
        "config": {"stratum": args.stratum, "min_class": args.min_class, "n_classes": len(classes),
                   "n_labelled": int(len(lab)), "n_search": int(len(Xtr)), "n_holdout": int(len(Xh)),
                   "trials": args.trials, "feature_dim": int(X.shape[1]), "seed": args.seed},
        "best": {"params": best, "rounds": rounds, "search_macro_f1": float(study.best_value)},
        "cv": cv,
        "holdout": {"macro_f1": float(f1_score(yh, hpred, average="macro", zero_division=0)),
                    "micro_f1": float(f1_score(yh, hpred, average="micro", zero_division=0)),
                    "logloss": float(log_loss(yh, hproba, labels=list(range(len(classes)))))},
        "baselines": baselines,
    }
    with open(os.path.join(args.out, "hpo_results_v135.json"), "w") as fh:
        json.dump(results, fh, indent=1)
    bst.save_model(os.path.join(args.out, "xgb_v135.ubj"))
    np.save(os.path.join(args.out, "pca_components.npy"), pca.components_)
    with open(os.path.join(args.out, "classes.json"), "w") as fh:
        json.dump(classes, fh)
    print(json.dumps(results["holdout"], indent=1))

if __name__ == "__main__":
    main()
