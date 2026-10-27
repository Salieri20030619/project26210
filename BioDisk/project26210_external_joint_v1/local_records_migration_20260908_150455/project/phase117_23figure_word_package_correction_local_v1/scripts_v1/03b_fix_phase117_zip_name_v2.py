from __future__ import annotations

import json
import shutil
from pathlib import Path


DOWNLOAD = Path(r"C:\Users\31937\Downloads\项目26210_2.1-2.2图件证据包_23图完整版_20260902_v2")
WRONG_ZIP = DOWNLOAD.with_suffix(".zip")
RIGHT_ZIP = Path(str(DOWNLOAD) + ".zip")
PHASE = Path(r"D:\project26210_external_joint_v1\phase117_23figure_word_package_correction_local_v1")
VALIDATION = PHASE / "INSTALL_AND_DELIVERY_VALIDATION_v2.json"

if not WRONG_ZIP.exists():
    raise FileNotFoundError(WRONG_ZIP)
if RIGHT_ZIP.exists():
    raise FileExistsError(RIGHT_ZIP)
WRONG_ZIP.replace(RIGHT_ZIP)

data = json.loads(VALIDATION.read_text(encoding="utf-8"))
data["zip"] = str(RIGHT_ZIP)
VALIDATION.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")

shutil.copy2(Path(__file__), PHASE / "scripts_v1" / Path(__file__).name)
print(str(RIGHT_ZIP))
