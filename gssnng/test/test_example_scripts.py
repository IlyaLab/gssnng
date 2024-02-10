
import time
import pathlib
import runpy

def run_examples():
    t0 = time.time()
    print('start time: ' + str(t0))
    # grabbing the example scripts
    scriptPath = pathlib.Path('.')
    scripts = list(scriptPath.glob('example_*.py'))
    
    print("RUNNING SCRIPTS: ")
    for si in scripts:
        print(si)
        runpy.run_path(si)
    print("")

    print('******DONE*******')
    t1 = time.time()
    print('end time: ' + str(t1))
    print('TOTAL TIME: ' + str(t1-t0))
    print('test done')
    return(True)

def test_that_examples_run():
    assert run_examples() == True  ## The up and dn sets should be combined into one.

