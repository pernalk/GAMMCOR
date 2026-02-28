from scripts.colors import *

def convert_level_to_int(level):
    levels = {
        'short': 0,
        'long': 1,
        'full': 2
    }
    return levels.get(level)

def should_run(level, TEST_LEVEL):
    level_int = convert_level_to_int(level)
    TEST_LEVEL_int = convert_level_to_int(TEST_LEVEL)
    if level_int <= TEST_LEVEL_int:
        return True
    else:
        return False

def prepare_test(test, TEST_LEVEL):
    units = []
    for unit, details in test.get('units').items():
        if should_run(details.get('level'), TEST_LEVEL):
            units.append({'name': unit,
                'atol': details.get('atol') })
    return units

def start_test(test, units, RUN_GAMMCOR):
    print(f"{BOLD}{test.get('name')}{ENDC}")
    runner = test.get('runner')
    return runner(units, test.get('fun'), RUN_GAMMCOR)