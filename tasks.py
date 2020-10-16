""" Some tasks for deployment to PyPI/TestPyPI. """


from invoke import task


@task
def publish_to_pypi(c, test=False):
    c.run("rm -rf build dist", warn=True)
    c.run("python setup.py sdist bdist_wheel")
    pypi_cmd = "twine upload --verbose dist/*"

    if test:
        pypi_cmd += " --repository-url https://test.pypi.org/legacy/"

    c.run(pypi_cmd, pty=True)
    c.run("rm -rf build dist", warn=True)
