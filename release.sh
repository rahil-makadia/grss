server=$1
if [ "$server" = "-pypi" ]; then
    python3 -m twine upload dist/*
else
    python3 -m twine upload --repository testpypi dist/*
fi

# install from pypi or testpypi based on -pypi flag
if [ "$server" = "-pypi" ]; then
    python3 -m pip install --upgrade grss
else
    python3 -m pip install --upgrade --index-url https://test.pypi.org/simple/ --no-deps grss --extra-index-url https://pypi.org/simple
fi
