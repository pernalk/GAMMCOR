#!/usr/bin/env python3

import argparse
import sys
import h5py


def display_value(value):
    if isinstance(value, bytes):
        value = value.decode(errors="replace")
    if isinstance(value, str) and (len(value) > 120 or "\n" in value):
        return f"<string, {len(value)} characters>"
    return value


def describe(name, obj):
    indent = "  " * name.count("/")
    if isinstance(obj, h5py.Dataset):
        print(f"{indent}[Dataset] {name}  shape={obj.shape}  dtype={obj.dtype}")
    elif isinstance(obj, h5py.Group):
        print(f"{indent}[Group]   {name}/")

    for attr_name, attr_val in obj.attrs.items():
        print(f"{indent}    @{attr_name} = {display_value(attr_val)}")


def print_selected(h5file, requested):
    name = requested.strip("/")
    if name not in h5file:
        print(f"Dataset not found: {requested}", file=sys.stderr)
        return False

    obj = h5file[name]
    if isinstance(obj, h5py.Group):
        describe(name, obj)
        obj.visititems(lambda child_name, child: describe(f"{name}/{child_name}", child))
        return True

    value = obj[()]
    if getattr(value, "size", 0) == 1:
        value = value.reshape(-1)[0]
        if hasattr(value, "item"):
            value = value.item()
    print(f"{name} = {display_value(value)}")
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("hdf5")
    parser.add_argument("dataset", nargs="*")
    args = parser.parse_args()

    with h5py.File(args.hdf5, "r") as f:
        if args.dataset:
            return 0 if all(print_selected(f, name) for name in args.dataset) else 1

        print(f"=== Structure of {args.hdf5} ===\n")

        for attr_name, attr_val in f.attrs.items():
            print(f"@{attr_name} = {display_value(attr_val)}")
        if f.attrs:
            print()

        f.visititems(describe)

    return 0


if __name__ == "__main__":
    sys.exit(main())
