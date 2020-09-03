If ($env:APPVEYOR_REPO_TAG -eq "true") {
    Invoke-Expression "twine upload --skip-existing dist/*.whl" 2>$null
} Else {
    write-output "Not on a tag on master, won't deploy to PyPI"
}
