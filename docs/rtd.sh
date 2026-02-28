. venv/bin/activate
python -m pip install --upgrade --no-cache-dir pip setuptools
python -m pip install --exists-action=w --no-cache-dir -r requirements.txt
python -m sphinx -T -b html -d _build/doctrees -D language=en . html
