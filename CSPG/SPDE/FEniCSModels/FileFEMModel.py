from FEMModel import *
parameters["refinement_algorithm"] = "plaza_with_parent_facets" # Bug fix for dolfin 1.5

import os


class FileFEMModel(FEMModel):
    def init_mesh(self, file_prefix):
        __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

        self.mesh       = Mesh(os.path.join(__location__, file_prefix + ".xml"))
        self.subdomains = MeshFunction("size_t", self.mesh, os.path.join(__location__, file_prefix + "_physical_region.xml"))
        self.boundaries = MeshFunction("size_t", self.mesh, os.path.join(__location__, file_prefix + "_facet_region.xml"))

    def __getstate__(self):
        odict = self.__dict__.copy()

        # Can't pickle mesh, solver and M
        del odict['mesh']
        del odict['subdomains']
        del odict['boundaries']
        del odict['solver']
        del odict['M']

        return odict
