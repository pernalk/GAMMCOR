import os
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

def get_interface(filename):
    with open(filename, 'r') as file:
        for line in file:
            words = line.lower().split()
            if len(words) > 1 and words[0] == "interface":
                return words[1]
    return None

def set_basis_path(filenames, basis_path):
    for filename in filenames:
        append_string_to_line(filename, "BasisPath", basis_path)

def clear_basis_path(filenames):
    for filename in filenames:
        remove_text_from_line(filename, "BasisPath")

def check_word_in_file(filename, word):
    with open(filename, 'r') as file:
        for line in file:
            if word.lower() in line.lower():
                return True
    return False

def remove_text_from_line(filename, line_prefix):
    with open(filename, 'r') as file:
        lines = file.readlines()

    with open(filename, 'w') as file:
        for line in lines:
            if line.lower().lstrip().startswith(line_prefix.lower()):
                line = line.split(maxsplit=1)[0]+'\n'
                line = " "+line
            file.write(line)

def append_string_to_line(filename, word, value):
    with open(filename, 'r') as file:
        lines = file.readlines()

    # append with value
    with open(filename, 'w') as file:
        for line in lines:
            if line.lower().lstrip().startswith(word.lower()):
                line = line.rstrip('\n') + ' ' + value + '\n'
            file.write(line)

