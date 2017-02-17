.. Documentation to the GraGLeS 2D Project
.. ===================================

The Grain Growth Level Set simulation tool is a 2D and 3D energy dissipation solver for arbitrary tesselations. The evolution is modeled as a competition of local gradient flows of each phase. 

The physical phenomena of grain growth in polycrystalline materials is the targed application. Thus the influence of mena curcvature flow as well as bulk energies is considered. The GraGLes project is developed by Christian Mießen and Co-Workers at the Institute of Physical Metallurgy and Metal Physics.

Please follow the installation instructions:

 .. figure:: ../images/IMMLogo.png
   :scale: 50%
   :align: center


1. Getting started
^^^^^^^^^^^^^^^^^^

The compilation of GraGLeS utilizes standard Linux tools.

.. toctree::
   :maxdepth: 2
     
   setup
   


2. Visualizing results in 2D
^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2
     
   visualize
   


3. Topological Evolution
^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2
     
   topotracer
   





References
^^^^^^^^^^

Specifically about GraGLeS:
 | Miessen C., Liesenjohann M., Barrales-Mora L.A., Shvindlerman L.S., Gottstein G.
 | An advanced level set approach to grain growth – Accounting for grain boundary anisotropy and finite triple junction mobility
 | Acta Materialia, 2015, 99
 | doi:10.1016/j.actamat.2015.07.040
 
 | Christian Mießen, Nikola Velinov, Günter Gottstein, Luis A. Barrales-Mora
 | A highly efficient 3D level-set grain growth algorithm tailored for ccNUMA architecture
 | https://arxiv.org/abs/1701.06658
 
 
For the post-analysis:
 | Kuhbach M., Barrales-Mora L.A., Miessen C., Gottstein G.: 
 | **Ultrafast analysis of individual grain behavior during grain growth by parallel computing** 
 | Proceedings of the 36th Riso International Symposium on Materials Science 
 | doi:10.1088/1757-899X/89/1/012031
 
Funding:
 | The authors gratefully acknowledge the support from the DFG in the frame of the Reinhart Kosselleck project (GO 335/44-1).
 | Furthermore, we acknowledge the support from the FZJuelich and RWTH Aachen University within the project JARAHPC projects. 
 
 
Version history
^^^^^^^^^^^^^^^

 | **v0.1** 

   
Licence
^^^^^^^

The project is licenced under the GNU v2.0


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

Questions, contributions
^^^^^^^^^^^^^^^^^^^^^^^^

Just let me know or contact *miessen@imm.rwth-aachen.de*
