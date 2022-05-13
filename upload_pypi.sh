# ====== build package ======
rm dist/*
rm *.egg-info/*
python3 -m build

# upload
pip install twine --upgrade
python3 -m twine upload --repository pypi dist/*

# change version info:
# 1. setup.py
# 2. __init__.py

# name: Cynthxxx
# pw: *****-Nbs

