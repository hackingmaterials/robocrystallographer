"""
Deployment file to facilitate releases of robocrystallographer.
Note that this file is meant to be run from the root directory of the repo.
"""

from invoke import task

import os
import json
import requests
import re
from robocrys import __version__
from monty.os import cd

exclude_paths = "'../robocrys/*/tests' ../robocrys/tests"


@task
def make_doc(ctx):
    with cd("docs_rst"):
        ctx.run("sphinx-apidoc -f -o source ../robocrys " + exclude_paths)
        ctx.run("make html")
        ctx.run("cp -r build/html/* ../docs")

    with cd("docs"):
        # Avoid the use of jekyll so that _dir works as intended.
        ctx.run("touch .nojekyll")


@task
def update_doc(ctx):
    make_doc(ctx)
    with cd("docs"):
        ctx.run("git add .")
        ctx.run("git commit -a -m \"Update to v{}\"".format(__version__))
        ctx.run("git push")

@task
def publish(ctx):
    ctx.run("rm dist/*.*", warn=True)
    ctx.run("python3 setup.py sdist bdist_wheel")
    ctx.run("twine upload dist/*")


@task
def release(ctx):
    with open("CHANGELOG.rst") as f:
        contents = f.read()
    new_ver = re.findall('\n(v.*)', contents)[0]
    toks = re.finditer("v\d\.\d\.\d*\n\-*(.*?)^v\d\.\d\.\d", contents,
                       re.MULTILINE | re.DOTALL)
    desc = list(toks)[0].groups()[0].strip()
    payload = {
        "tag_name": new_ver,
        "target_commitish": "master",
        "name": new_ver,
        "body": desc,
        "draft": False,
        "prerelease": False
    }
    response = requests.post(
        "https://api.github.com/repos/hackingmaterials/robocrystallographer/releases",
        data=json.dumps(payload),
        headers={"Authorization": "token " + os.environ["GITHUB_TOKEN"]})
    print(response.text)
