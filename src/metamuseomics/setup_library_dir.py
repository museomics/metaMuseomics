```python
import shutil
import subprocess
from pathlib import Path
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed


def run_cmd(cmd, cwd=None):
    """Run a command and return (success: bool)."""
    try:
        subprocess.run(cmd, cwd=cwd, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {' '.join(cmd)} in {cwd}\n{e}")
        return False


def download_file(url: str, output_dir: Path):
    """Download a file with resume support."""
    print(f"Downloading: {url}")
    try:
        subprocess.run(
            ["wget", "-c", url, "-P", str(output_dir)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        print(f"Finished: {url}")
    except subprocess.CalledProcessError as e:
        print(f"Download failed for {url}: {e}")
        raise


def main():
    parser = argparse.ArgumentParser(
        description="Recoverable downloader / extractor."
    )
    parser.add_argument("master_dir", type=str)
    parser.add_argument("links_file", type=str)
    parser.add_argument("--url_id", type=str, default=None)
    parser.add_argument("--threads", type=int, default=4)

    args = parser.parse_args()

    master_dir = Path(args.master_dir)
    links_file = Path(args.links_file)

    # READ URL LIST (OPTIONAL)
    if links_file.exists() and links_file.stat().st_size > 0:
        with links_file.open() as f:
            urls = [x.strip() for x in f if x.strip()]
        print(f"Found {len(urls)} URLs.")
    else:
        print("Links file empty or missing — skipping download phase.")
        urls = []

    # Determine directory ID
    url_id = args.url_id or "download"

    sub_dir = master_dir / url_id
    raw_data_dir = sub_dir / "raw_data"
    raw_data_dir.mkdir(parents=True, exist_ok=True)

    # 1. DOWNLOAD (optional)
    if urls:
        print(f"Downloading using {args.threads} threads…")

        with ThreadPoolExecutor(max_workers=args.threads) as ex:
            futures = {
                ex.submit(download_file, url, sub_dir): url
                for url in urls
            }

            for future in as_completed(futures):
                url = futures[future]
                try:
                    future.result()
                except Exception:
                    print(f"Failed downloading: {url}")
    else:
        print("No URLs provided — skipping download.")

    # 2. TOP-LEVEL MD5 CHECK
    md5_top = sub_dir / "MD5.txt"

    if md5_top.exists():
        print("\nRunning top-level MD5 validation (for tar files)…")

        if not run_cmd(
            ["md5sum", "-c", "MD5.txt"],
            cwd=sub_dir,
        ):
            print(
                "⚠ Some top-level MD5 checks failed. "
                "You may need to re-download or verify manually."
            )
    else:
        print("No top-level MD5.txt found — skipping MD5 verification.")

    # 3. EXTRACT TAR FILES
    print("\nExtracting tar files…")

    for tar_file in sub_dir.glob("*.tar"):
        extract_dir = sub_dir / tar_file.stem

        if extract_dir.exists():
            print(
                f"Skipping extraction for {tar_file} "
                "(directory already exists)."
            )
            continue

        print(f"Extracting {tar_file}…")

        if run_cmd(
            [
                "tar",
                "-xvf",
                str(tar_file),
                "-C",
                str(sub_dir),
            ]
        ):
            tar_file.unlink()
            print(f"Removed {tar_file}.")
        else:
            print(
                f"Extraction failed for {tar_file}. Continuing."
            )

    # 4. INNER MD5 VALIDATION
    print("\nRunning inner MD5 checks…")

    for extracted_dir in sub_dir.iterdir():
        if (
            not extracted_dir.is_dir()
            or extracted_dir.name == "raw_data"
        ):
            continue

        inner_md5 = extracted_dir / "MD5.txt"

        if inner_md5.exists():
            print(
                f"Validating MD5 inside {extracted_dir}…"
            )
            run_cmd(
                ["md5sum", "-c", "MD5.txt"],
                cwd=extracted_dir,
            )
        else:
            print(
                f"No inner MD5 found inside {extracted_dir} "
                "— skipping."
            )

    # 5. MOVE fq.gz FILES INTO raw_data/
    print("\nCollecting .fq.gz files…")

    for fq in sub_dir.rglob("*.fq.gz"):
        dst = raw_data_dir / fq.name

        if fq != dst:
            fq.rename(dst)
            print(f"Moved {fq} → {dst}")

    # 6. CLEANUP — remove leftover extracted folders
    print("\nCleaning up extracted directories…")

    for p in sub_dir.iterdir():
        if p.is_dir() and p != raw_data_dir:
            shutil.rmtree(p)
            print(f"Removed {p}")

    print("\nCOMPLETE")
    print(f"Final processed files stored in:\n{raw_data_dir}")


if __name__ == "__main__":
    main()
```
