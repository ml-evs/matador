""" Some tasks for deployment to PyPI/TestPyPI. """


from invoke import task


@task
def publish_to_pypi(c, test=False):
    c.run("rm -rf build dist", warn=True)
    c.run("python setup.py sdist bdist_wheel")
    c.run("twine upload --verbose dist/*" + " --repository-url https://test.pypi.org/legacy/" if test else "")
    c.run("rm -rf build dist", warn=True)
