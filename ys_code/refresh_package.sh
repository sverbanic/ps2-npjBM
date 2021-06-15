#!/bin/bash

python setup.py bdist_wheel
pip uninstall skin_mb_chen
MB_PKG=$(ls -t dist/*.whl | head -1) 
pip install $MB_PKG
