import re
import setuptools

# function to read the version wrote down in the project_name.__init__ file
def get_property(prop, project):
    result = re.search(r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop), open(project + '/__init__.py').read())
    return result.group(1)


project_name = 'aramis_exporter'

setuptools.setup(
    name='aramis_exporter',
    version=get_property('__version__', project_name),
    packages=setuptools.find_packages(
        exclude=['test_scripts*',
                 'example_images']
    ),
    description='Zeiss GOM Aramis Data Exporter',
    author='DLR',
    license='MIT',
    include_package_data=True,
    install_requires=[
        'numpy',
    ]
)
