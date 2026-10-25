"""Wrap a relocatable runtime archive in a self-extracting installer."""

import argparse
import hashlib
import shutil
from pathlib import Path

HEADER = r"""#!/usr/bin/env bash
set -euo pipefail
prefix="$HOME/.local/share/openaschem/1.3.0"
while (($#)); do
  case "$1" in
    --prefix) [[ $# -ge 2 ]] || { echo "--prefix requires a directory"; exit 2; }; prefix="$2"; shift 2 ;;
    -h|--help) echo 'Usage: bash OpenASChem-1.3.0-Linux-x86_64.sh [--prefix DIRECTORY]'; echo 'Includes Python and all dependencies. No downloads required.'; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done
[[ "$(uname -m)" == x86_64 ]] || { echo 'This installer requires Linux x86_64.' >&2; exit 2; }
for utility in tar sha256sum tail awk tee; do
  command -v "$utility" >/dev/null || { echo "Missing required command: $utility" >&2; exit 1; }
done
mkdir -p -- "$prefix"
prefix=$(cd -- "$prefix" && pwd)
if [[ -d "$prefix/environment" ]]; then
  [[ -f "$prefix/openaschem-install.json" ]] || { echo 'Choose a new empty installation directory.' >&2; exit 1; }
  metadata=$(<"$prefix/openaschem-install.json")
  [[ "$metadata" == *bundled* || "$metadata" == *fully-offline* ]] || { echo 'Choose a new installation directory.' >&2; exit 1; }
fi
archive_line=$(awk '/^__OPENASCHEM_ENVIRONMENT__$/ {print NR+1; exit}' "$0")
echo 'Checking bundled environment...'
actual=$(tail -n +"$archive_line" "$0" | sha256sum)
[[ "${actual%% *}" == '__SHA256__' ]] || { echo 'Installer checksum mismatch.' >&2; exit 1; }
printf '%s\n' '{"application":"OpenASChem","version":"1.3.0","status":"installing","installation_mode":"bundled"}' > "$prefix/openaschem-install.json"
mkdir -p -- "$prefix/environment"
echo 'Extracting Python and dependencies. No downloads are required.'
tail -n +"$archive_line" "$0" | tar -xzf - -C "$prefix/environment"
"$prefix/environment/bin/python" -E "$prefix/environment/openaschem-finalize.py" --prefix "$prefix" 2>&1 | tee "$prefix/installation.log"
echo "GUI: $prefix/bin/openaschem-gui"
echo "CLI: $prefix/bin/openaschem"
exit 0
__OPENASCHEM_ENVIRONMENT__
"""


def build(archive, output):
    checksum = hashlib.sha256(archive.read_bytes()).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as target:
        target.write(HEADER.replace("__SHA256__", checksum).encode())
        with archive.open("rb") as stream:
            shutil.copyfileobj(stream, target)
    output.chmod(0o755)
    print(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument(
        "--out", type=Path, default=Path("release/OpenASChem-1.3.0-Linux-x86_64.sh")
    )
    args = parser.parse_args()
    build(args.archive, args.out)
