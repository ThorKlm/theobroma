"""Download THEOBROMA dataset from HuggingFace."""
import argparse, os

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--token", default=os.environ.get("HF_TOKEN", ""))
    ap.add_argument("--out", default="/home/thorben.klamt/theobroma/data")
    ap.add_argument("--repo", default="ThorKl/theobroma")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    from huggingface_hub import snapshot_download
    snapshot_download(repo_id=args.repo, repo_type="dataset",
                      local_dir=args.out, token=args.token or None)
    print(f"Dataset downloaded to {args.out}")
    for f in sorted(os.listdir(args.out)):
        size = os.path.getsize(os.path.join(args.out, f)) / (1024*1024)
        print(f"  {f} ({size:.1f} MB)")

if __name__ == "__main__":
    main()
