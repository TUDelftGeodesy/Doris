from setuptools import setup

setup(
    name='doris',
    version='5.0.3',
    packages=[
        'doris',
        'doris.doris_stack',
        'doris.doris_stack.functions',
        'doris.doris_stack.main_code',
        'doris.prepare_stack',
    ],
    package_dir = {
        'doris': '',
        'doris.doris_stack': 'doris_stack',
        'doris.doris_stack.functions': 'doris_stack/functions',
        'doris.doris_stack.main_code': 'doris_stack/main_code',
        'doris.prepare_stack': 'prepare_stack',
    },
    url='https://github.com/TUDelftGeodesy/Doris',
    license='LICENSE.txt',
    author='Gert Mulder',
    author_email='g.mulder-@tudelft.nl',
    description='doris InSAR processing software',
    install_requires=['numpy', 'shapely', 'requests', 'fiona', 'gdal', 'osr', 'scipy', 'fastkml']
)
