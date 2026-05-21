#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
env_file="$script_dir/environment-python-spacel3d.yml"

if conda env list | awk '{print $1}' | grep -qx 'stereo-skills-py-spacel3d'; then
  echo "Conda env stereo-skills-py-spacel3d already exists; updating base dependencies."
  conda env update -n stereo-skills-py-spacel3d -f "$env_file" --prune
else
  conda env create -f "$env_file"
fi

conda run -n stereo-skills-py-spacel3d python -m pip install --no-deps SPACEL==1.1.8
conda run -n stereo-skills-py-spacel3d python -c "import importlib.util; from pathlib import Path; import SPACEL; plot_path = Path(SPACEL.__file__).resolve().parent / 'Scube' / 'plot.py'; spec = importlib.util.spec_from_file_location('spacel_scube_plot', plot_path); module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module); assert hasattr(module, 'plot_3d'); print(f'SPACEL Scube plotting import OK: {plot_path}')"
