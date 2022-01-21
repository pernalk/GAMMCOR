#!/bin/bash

TREXIO_PATH=$1
LIBS=("./libtrexio.so" "./libtrexio.so.*" "./libtrexio_f.so" "../SOURCE/trexio_f.f90")

echo "Copying files..."

cp -f ${TREXIO_PATH}/build/src/libtrexio*.so* .
cp -f ${TREXIO_PATH}/include/trexio_f.f90 ../SOURCE/

echo "Verifying..."

status=0
for file in ${LIBS[@]}; do
  if test -f ${file}; then
    echo -e "\e[32m[x] ${file}\e[0m"
  else
    echo -e "\e[31m[ ] ${file}\e[0m"
    status=-1
  fi
done

if [ ${status} == 0 ]; then
  echo -e "\e[32mInstallation complete!\e[0m"
else
  echo -e "\e[31mInstallation failed!\e[0m"
fi
