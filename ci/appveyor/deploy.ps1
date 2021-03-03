If ($env:APPVEYOR_REPO_TAG -eq "true") {
    Invoke-Expression "$env:PYTHON\\python.exe -m pip install twine"
    Invoke-Expression "$env:PYTHON\\python.exe -m twine upload --skip-existing dist\\*.whl"
} Else {
    write-output "Not on a tag on master, won't deploy to PyPI"
}
