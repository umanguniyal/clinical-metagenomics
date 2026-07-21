#!/usr/bin/env python3
import argparse, os, re, subprocess, sys, time
from pathlib import Path
from urllib.parse import urlparse
import gzip

ENA_HTTP = "https://ftp.sra.ebi.ac.uk/vol1/fastq"

def is_url(s: str) -> bool:
    return s.startswith(("http://", "https://", "ftp://"))

def is_accession(s: str) -> bool:
    return bool(re.match(r"^(ERR|SRR|DRR|ERX|SRX|DRX)\d+$", s.strip()))

def ena_base(acc: str) -> str:
    prefix = acc[:6]
    subdir = ""
    if len(acc) > 9:
        try:
            subdir = "/" + f"{int(acc[9:]):03d}"
        except:
            subdir = ""
    return f"{ENA_HTTP}/{prefix}{subdir}/{acc}"

def aria2c_download(url: str, outdir: str, max_retries: int = 3) -> bool:
    os.makedirs(outdir, exist_ok=True)
    for attempt in range(max_retries):
        try:
            print(f"Aria2c attempt {attempt + 1}/{max_retries}")
            cmd = [
                "aria2c",
                "-x", "16",
                "-s", "16",
                "-j", "4",
                "--retry-wait=2",
                "--max-tries=5",
                "--timeout=120",
                "--check-integrity=true",
                "-d", outdir,
                url
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            if result.returncode == 0:
                print(f"✓ Aria2c successful")
                return True
            else:
                if attempt < max_retries - 1:
                    time.sleep(3)
        except Exception as e:
            if attempt < max_retries - 1:
                time.sleep(3)
    return False

def curl_download(url: str, outpath: str, max_retries: int = 3) -> bool:
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    for attempt in range(max_retries):
        try:
            print(f"Curl attempt {attempt + 1}/{max_retries}")
            cmd = [
                "curl",
                "-L",
                "-C", "-",
                "--connect-timeout", "30",
                "--max-time", "300",
                "--retry", "3",
                "--retry-delay", "2",
                "-o", outpath,
                url
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            if result.returncode == 0 and os.path.exists(outpath) and os.path.getsize(outpath) > 1000000:
                print(f"✓ Curl successful")
                return True
            else:
                if attempt < max_retries - 1:
                    time.sleep(3)
        except Exception as e:
            if attempt < max_retries - 1:
                time.sleep(3)
    return False

def validate_fastq_file(path: str) -> bool:
    """Check gzip integrity and FASTQ @ record."""
    try:
        if path.endswith(".gz"):
            f = gzip.open(path, "rt", errors="ignore")
        else:
            f = open(path, "rt", errors="ignore")
        found = False
        n = 0
        for line in f:
            if line.startswith("@"):
                found = True
                break
            n += 1
            if n > 5:
                break
        f.close()
        if not found:
            print(f"FAIL: {path} - no '@' found in first record", file=sys.stderr)
            return False
        if os.path.getsize(path) < 10000:
            print(f"FAIL: {path} - file is very small", file=sys.stderr)
            return False
    except Exception as e:
        print(f"FAIL: {path} - cannot read/parse ({e})", file=sys.stderr)
        return False
    return True

def download_file(url: str, outpath: str) -> bool:
    """
    Download a single file, trying aria2c first, then curl, then validate FASTQ.
    """
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    if aria2c_download(url, os.path.dirname(outpath)):
        if not validate_fastq_file(outpath):
            print(f"✗ ERROR: FASTQ failed validation: {outpath}", file=sys.stderr)
            if os.path.exists(outpath):
                os.remove(outpath)
            return False
        return True
    print(f"Aria2c failed, trying curl...")
    if curl_download(url, outpath):
        if not validate_fastq_file(outpath):
            print(f"✗ ERROR: FASTQ failed validation: {outpath}", file=sys.stderr)
            if os.path.exists(outpath):
                os.remove(outpath)
            return False
        return True
    return False

def normalize_fastq_gz(path: str) -> str:
    p = Path(path)
    if p.name.endswith(".fastq.gz") or p.name.endswith(".fq.gz"):
        return str(p)
    if p.name.endswith(".fastq") or p.name.endswith(".fq"):
        gz = str(p) + ".gz"
        if not os.path.exists(gz):
            subprocess.check_call(["bash", "-lc", f"gzip -c '{p}' > '{gz}'"])
        return gz
    return str(p)

def stage_symlink(src: str, dst: str):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.exists(dst):
        return
    os.symlink(os.path.abspath(src), dst)

def parse_fastq_name(name: str):
    m = re.search(r"_R?([12])(?:_\d+)?\.(?:fastq|fq)(?:\.gz)?$", name)
    if m:
        base = name[:m.start()]
        base = re.sub(r"_S\d+_L\d+$", "", base)
        return base, int(m.group(1))
    
    m2 = re.search(r"\.(?:fastq|fq)(?:\.gz)?$", name)
    if m2:
        return name[:m2.start()], None
        
    return None, None

def url_basename(url: str) -> str:
    return os.path.basename(urlparse(url).path)

def download_accession_illumina(acc: str, project_root: str):
    """
    Download PAIRED-END (Illumina) sample
    Must have both R1 and R2
    """
    base = ena_base(acc)
    illum_dir = os.path.join(project_root, "illumina", acc)
    os.makedirs(illum_dir, exist_ok=True)
    r1_url = f"{base}/{acc}_1.fastq.gz"
    r2_url = f"{base}/{acc}_2.fastq.gz"
    r1_path = os.path.join(illum_dir, f"{acc}_1.fastq.gz")
    r2_path = os.path.join(illum_dir, f"{acc}_2.fastq.gz")
    print(f"\n{'='*70}")
    print(f"ILLUMINA PAIRED-END: {acc}")
    print(f"{'='*70}")
    print(f"\nDownloading R1...")
    if not download_file(r1_url, r1_path):
        print(f"✗ ERROR: Failed to download R1")
        return None
    print(f"✓ R1 downloaded")
    print(f"\nDownloading R2...")
    if not download_file(r2_url, r2_path):
        print(f"✗ ERROR: Failed to download R2")
        print(f"Note: This is nanopore data (single-end). Use platform 2 (nanopore)")
        return None
    print(f"✓ R2 downloaded")
    print(f"\n✓ SUCCESS: {acc} - ILLUMINA (paired-end)")
    return ("illumina", acc, r1_path, r2_path)

def download_accession_nanopore(acc: str, project_root: str):
    """
    Download SINGLE-END (Nanopore) sample
    Only tries R1 (or .fastq without number suffix)
    """
    base = ena_base(acc)
    nano_dir = os.path.join(project_root, "nanopore", acc)
    os.makedirs(nano_dir, exist_ok=True)
    r1_url = f"{base}/{acc}_1.fastq.gz"
    r1_path = os.path.join(nano_dir, f"{acc}_1.fastq.gz")
    print(f"\n{'='*70}")
    print(f"NANOPORE SINGLE-END: {acc}")
    print(f"{'='*70}")
    print(f"\nDownloading single-end read...")
    if download_file(r1_url, r1_path):
        print(f"✓ Single-end read downloaded")
        print(f"\n✓ SUCCESS: {acc} - NANOPORE (single-end)")
        return ("nanopore", acc, r1_path, "")
    alt_url = f"{base}/{acc}.fastq.gz"
    alt_path = os.path.join(nano_dir, f"{acc}.fastq.gz")
    print(f"\nTrying alternative naming: {url_basename(alt_url)}")
    if download_file(alt_url, alt_path):
        os.rename(alt_path, r1_path)
        print(f"✓ Read downloaded and renamed")
        print(f"\n✓ SUCCESS: {acc} - NANOPORE (single-end)")
        return ("nanopore", acc, r1_path, "")
    print(f"✗ ERROR: Failed to download single-end read for {acc}")
    return None

def scan_staged(project_root: str):
    candidates = {}
    def add_file(platform, fpath):
        bn = os.path.basename(fpath)
        sid, mate = parse_fastq_name(bn)
        if not sid:
            return
        if sid not in candidates:
            candidates[sid] = {"illumina": {}, "nanopore": {}}
        if platform == "illumina":
            candidates[sid]["illumina"][mate] = fpath
        else:
            candidates[sid]["nanopore"][mate] = fpath
    illum_root = os.path.join(project_root, "illumina")
    nano_root = os.path.join(project_root, "nanopore")
    if os.path.isdir(illum_root):
        for sample_dir in os.listdir(illum_root):
            sample_path = os.path.join(illum_root, sample_dir)
            if not os.path.isdir(sample_path):
                continue
            for ext in ("*.fastq*", "*.fq*"):
                for fpath in Path(sample_path).rglob(ext):
                    if "trimmed" in fpath.name or "nonhost" in fpath.name:
                        continue
                    add_file("illumina", str(fpath))
    if os.path.isdir(nano_root):
        for sample_dir in os.listdir(nano_root):
            sample_path = os.path.join(nano_root, sample_dir)
            if not os.path.isdir(sample_path):
                continue
            for ext in ("*.fastq*", "*.fq*"):
                for fpath in Path(sample_path).rglob(ext):
                    if "trimmed" in fpath.name or "nonhost" in fpath.name:
                        continue
                    add_file("nanopore", str(fpath))
    rows = []
    for sid, info in candidates.items():
        illum_files = info["illumina"]
        nano_files = info["nanopore"]
        if 1 in illum_files and 2 in illum_files:
            r1 = normalize_fastq_gz(illum_files[1])
            r2 = normalize_fastq_gz(illum_files[2])
            outdir = os.path.join(project_root, "illumina", sid)
            stage_symlink(r1, os.path.join(outdir, f"{sid}_1.fastq.gz"))
            stage_symlink(r2, os.path.join(outdir, f"{sid}_2.fastq.gz"))
            rows.append(("illumina", sid,
                         os.path.join(outdir, f"{sid}_1.fastq.gz"),
                         os.path.join(outdir, f"{sid}_2.fastq.gz")))
            print(f"Sample {sid}: ILLUMINA (paired-end)")
        elif 1 in nano_files or None in nano_files:
            src = nano_files.get(1) or nano_files.get(None)
            if src:
                fq = normalize_fastq_gz(src)
                outdir = os.path.join(project_root, "nanopore", sid)
                stage_symlink(fq, os.path.join(outdir, f"{sid}_1.fastq.gz"))
                rows.append(("nanopore", sid,
                             os.path.join(outdir, f"{sid}_1.fastq.gz"),
                             ""))
                print(f"Sample {sid}: NANOPORE (single-end)")
        elif 1 in illum_files and 2 not in illum_files:
            src = illum_files[1]
            fq = normalize_fastq_gz(src)
            outdir = os.path.join(project_root, "nanopore", sid)
            os.makedirs(outdir, exist_ok=True)
            dst = os.path.join(outdir, f"{sid}_1.fastq.gz")
            stage_symlink(fq, dst)
            rows.append(("nanopore", sid, dst, ""))
            print(f"Sample {sid}: NANOPORE (single-end)")
    rows.sort(key=lambda x: (x[0], x[1]))
    return rows

def ingest_inputs(inputs, project_root: str, platform: str):
    candidates = {}
    def add_path(fpath, override_sid=None):
        bn = os.path.basename(fpath)
        sid, mate = parse_fastq_name(bn)
        if override_sid:
            sid = override_sid
        elif not sid:
            sid = re.sub(r"\.(?:fastq|fq)(?:\.gz)?$", "", bn)
            mate = None
        if not sid or sid.lower() in ("1", "2", "r1", "r2", "read", "reads"):
            sid = os.path.basename(os.path.dirname(fpath))
            
        if sid not in candidates:
            candidates[sid] = {"files": {}}
        if mate in (1, 2):
            candidates[sid]["files"][mate] = fpath
        else:
            candidates[sid]["files"][None] = fpath
    for item in inputs:
        item = item.strip()
        if not item:
            continue
        if is_accession(item):
            if platform == "illumina":
                result = download_accession_illumina(item, project_root)
            else:
                result = download_accession_nanopore(item, project_root)
            if result:
                plat, sid, r1, r2 = result
                if sid not in candidates:
                    candidates[sid] = {"files": {}}
                if plat == "illumina":
                    candidates[sid]["files"][1] = r1
                    candidates[sid]["files"][2] = r2
                else:
                    candidates[sid]["files"][1] = r1
            continue
        if is_url(item):
            tmp = os.path.join(project_root, "staging_downloads")
            os.makedirs(tmp, exist_ok=True)
            outfile = os.path.join(tmp, url_basename(item))
            if download_file(item, outfile):
                add_path(outfile)
            continue
        if os.path.isdir(item):
            for ext in ("*.fastq*", "*.fq*"):
                for p in Path(item).rglob(ext):
                    if "trimmed" in p.name or "nonhost" in p.name:
                        continue
                    add_path(str(p))
            continue
        if os.path.isfile(item):
            if re.search(r"\.(?:fastq|fq)(?:\.gz)?$", os.path.basename(item)):
                add_path(item)
            else:
                with open(item) as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith("#"):
                            continue
                        inputs.append(line)
            continue
        print(f"WARNING: Unrecognized input (skipped): {item}", file=sys.stderr)
    return finalize_candidates(candidates, project_root, platform)

def finalize_candidates(candidates, project_root: str, platform: str):
    rows = []
    for sid, info in candidates.items():
        files = info["files"]
        if 1 in files and 2 in files and platform == "illumina":
            r1 = normalize_fastq_gz(files[1])
            r2 = normalize_fastq_gz(files[2])
            outdir = os.path.join(project_root, "illumina", sid)
            stage_symlink(r1, os.path.join(outdir, f"{sid}_1.fastq.gz"))
            stage_symlink(r2, os.path.join(outdir, f"{sid}_2.fastq.gz"))
            rows.append(("illumina", sid,
                         os.path.join(outdir, f"{sid}_1.fastq.gz"),
                         os.path.join(outdir, f"{sid}_2.fastq.gz")))
        else:
            src = files.get(1) or files.get(None)
            if src:
                fq = normalize_fastq_gz(src)
                outdir = os.path.join(project_root, "nanopore", sid)
                stage_symlink(fq, os.path.join(outdir, f"{sid}_1.fastq.gz"))
                rows.append(("nanopore", sid,
                             os.path.join(outdir, f"{sid}_1.fastq.gz"),
                             ""))
    rows.sort(key=lambda x: (x[0], x[1]))
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    ap.add_argument("--out-manifest", required=True)
    ap.add_argument("--platform", required=True, choices=["illumina", "nanopore"])
    ap.add_argument("--inputs", nargs="*", default=[])
    args = ap.parse_args()
    project_root = os.path.abspath(args.project_root)
    os.makedirs(os.path.join(project_root, "illumina"), exist_ok=True)
    os.makedirs(os.path.join(project_root, "nanopore"), exist_ok=True)
    os.makedirs(os.path.dirname(args.out_manifest), exist_ok=True)
    if args.inputs and len(args.inputs) > 0:
        rows = ingest_inputs(list(args.inputs), project_root, args.platform)
    else:
        rows = scan_staged(project_root)
    with open(args.out_manifest, "w") as out:
        out.write("platform\tsample\tr1\tr2\n")
        for platform, sid, r1, r2 in rows:
            out.write(f"{platform}\t{sid}\t{r1}\t{r2}\n")
    print(f"\n{'='*70}")
    print(f"Manifest: {args.out_manifest}")
    print(f"Total samples: {len(rows)}")
    for platform, sid, r1, r2 in rows:
        status = "paired" if r2 else "single"
        print(f"  {platform}: {sid} ({status})")
    print(f"{'='*70}\n")

if __name__ == "__main__":
    main()
