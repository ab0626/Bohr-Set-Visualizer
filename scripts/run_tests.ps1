# Avoid conflicts from globally installed pytest plugins (e.g. web3).
$env:PYTEST_DISABLE_PLUGIN_AUTOLOAD = "1"
Set-Location $PSScriptRoot\..
python -m pytest tests/ -v --tb=short
