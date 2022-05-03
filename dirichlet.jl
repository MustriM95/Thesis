# Drawing a Dirichlet distribution inside a 2-simplex
using Random, Distributions, Plots, Meshes, MeshViz, GLMakie
p1 = Meshes.Point(0.0, 0.0)
p2 = Meshes.Point(1.0, 0.0)
p3 = Meshes.Point(0.5, 0.375)
p4 = Meshes.Point(-0.5, 0.375)
simp = Ngon( p1, p2, p3 )
quad = Ngon(p1, p2, p3, p4)
ref1 = refine(simp, TriRefinement())
ref2 = refine(ref1, CatmullClark())
ref3 = refine(ref2, CatmullClark())
ref4 = refine(ref3, TriRefinement())
viz(simp, showfacets = true)
viz(ref3, showfacets = true)
