# installer script
pip uninstall -y dsgrn_net_query &> /dev/null || True
python install_requirements.py
pip install -e .
