Python memory-efficient strings lookup
==============


Installation
------------

 - clone this repository
 - `python3 setup.py build`


To release
-----------

To install with dependencies:

```python3 -m pip install . ```

To generate a source distribution:

```python3 setup.py sdist ```

To upload the artifacts to the repo

```twine upload -r testpypi dist/* ```

To install from repo
--------

Run the following command:

```pip3 install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple cpp-string-low-mem-lookup==0.0.5```

License
-------

the package is provided under a BSD-style license that can be found in the LICENSE
file. By using, distributing, or contributing to this project, you agree to the
terms and conditions of this license.
