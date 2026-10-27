#!/usr/bin/env python3
import inspect
import json
from cell2location.models import RegressionModel, Cell2location

print(json.dumps({
    "RegressionModel.setup_anndata": str(inspect.signature(RegressionModel.setup_anndata)),
    "RegressionModel.__init__": str(inspect.signature(RegressionModel.__init__)),
    "RegressionModel.train": str(inspect.signature(RegressionModel.train)),
    "RegressionModel.export_posterior": str(inspect.signature(RegressionModel.export_posterior)),
    "Cell2location.setup_anndata": str(inspect.signature(Cell2location.setup_anndata)),
    "Cell2location.__init__": str(inspect.signature(Cell2location.__init__)),
    "Cell2location.train": str(inspect.signature(Cell2location.train)),
    "Cell2location.export_posterior": str(inspect.signature(Cell2location.export_posterior)),
}, indent=2))
