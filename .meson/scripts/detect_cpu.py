#!/usr/bin/env python3
import sys


def main():
    vendor = ''
    name = ''
    flags = set()
    with open('/proc/cpuinfo') as f:
        for line in f:
            key, _, value = line.partition(':')
            key = key.strip()
            if key == 'vendor_id' and not vendor:
                vendor = value.strip()
            elif key == 'model name' and not name:
                name = value.strip()
            elif key == 'flags' and not flags:
                flags = set(value.split())

    if {'avx512f', 'avx512cd', 'avx512bw', 'avx512dq', 'avx512vl'} <= flags:
        target = 'avx512' if vendor == 'GenuineIntel' else 'avx512-amd'
    elif {'avx2', 'fma'} <= flags:
        target = 'avx2'
    else:
        sys.exit('detect_cpu.py: this CPU has no AVX2/FMA, which is not supported')

    print(target)
    print(name)


if __name__ == '__main__':
    main()
