import os
import pytest
import oftest
import pandas as pd
from oftest import run_case
from oftest import run_reset_case
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

@pytest.fixture  # descruct all tests
def getProfilingData():
    profData = pd.DataFrame()
    dir_name = os.path.dirname(os.path.abspath(__file__))

    proc_dirs = [ f.path for f in os.scandir(dir_name) if "processor" in f.name]
    for proc in proc_dirs:
        d = ParsedParameterFile(f"{proc}/0/uniform/profiling")
        proc_data = pd.DataFrame(d["profiling"]).T
        proc_data["processor"] = proc.split("/")[-1]
        profData = profData.append(proc_data)
    return profData


def simMod(res):
    dir_name = os.path.dirname(os.path.abspath(__file__))
    filemod = {
        "system/simulationParameter": [("nx", res)]
    }
    meta_data = {"Res": res}
    case_mod = oftest.Case_modifiers(filemod, dir_name, meta_data)
    return case_mod


res25 = simMod(25)
res50 = simMod(50)
res100 = simMod(100)
res200 = simMod(200)



parameters = [
    res25,
    res50,
    res100,
    pytest.param(res200, marks=pytest.mark.slow)
]


@pytest.mark.parametrize(
    "run_reset_case",
    parameters,
    indirect=["run_reset_case"],
    ids=[
        "25",
        "50",
        "100",
        "200"
    ]
)
def test_buildstencil(run_reset_case, getProfilingData):
    log = oftest.path_log()
    dir_name = os.path.dirname(os.path.abspath(__file__))
    assert oftest.case_status(log) == "completed"  # checks if run completes

    grid_res = run_reset_case.meta_data["Res"]
    profData = getProfilingData
    masks = profData['description'].str.contains('calcCellStencil')
    profData["Res"] = grid_res
    profData = profData[masks]
    
    profData.to_csv(f"{dir_name}/profData{grid_res}.csv")


def test_write_results():
    dir_name = os.path.dirname(os.path.abspath(__file__))
    files = [ f.path for f in os.scandir(dir_name) if "profData" in f.name]
    profResults = pd.DataFrame()
    for f in files:
        profResults = profResults.append(pd.read_csv(f))

    profResults.to_csv(os.path.join(dir_name, "profDataTotal.csv"),index=False)