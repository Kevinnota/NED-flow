.. _instalation-page:

Instaling NED-flow
===========================

To use **NED-flow**, download the repository from GitHub:

.. code-block:: bash

   git clone https://github.com/kevinnota/NED-flow.git
   cd NED-flow


Python Dependencies
--------------------------

To use the Python tools in **NED-flow**, you must install the required Python packages. These are listed in `ned-py-install.txt`.

You can install them all at once using `pip`:

.. code-block:: bash

   pip install -r ned-py-install.txt

Alternatively, you can install the packages one by one:

.. code-block:: bash

   pip install tqdm==4.66.1
   pip install ete3==3.1.3
   pip install colorama==0.4.6
   pip install requests==2.31.0
   pip install pysam==0.22.0
   pip install scipy==1.12.0
   pip install biopython==1.83
   pip install numpy==1.26.3
   pip install intervaltree==3.1.0
   pip install scikit-learn==1.5.1
