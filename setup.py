from setuptools import setup, find_packages
from os import path
import versioneer

cur_dir = path.abspath(path.dirname(__file__))

# parse requirements
with open(path.join(cur_dir, "requirements.txt"), "r") as f:
    requirements = f.read().split()

setup(
    name="zfel",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    license='Apache Software License (http://www.apache.org/licenses/LICENSE-2.0)',
    license_files = ('LICENSE.txt',),
    install_requires=requirements,
    url="https://github.com/slaclab/zfel",
    include_package_data=True,
    python_requires=">=3.8",
)

