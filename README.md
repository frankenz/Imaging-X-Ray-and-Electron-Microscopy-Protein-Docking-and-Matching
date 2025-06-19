# Protein-protein docking and matching


Analyze the structural geometry, functional organization, dynamics and heterogeneity of Proteins. The procedures evaluate the geometric complementarity between the objects using grid representation and fast Fourier transformations.

Algorithm I developed published at:

https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0003594

https://www.cell.com/structure/fulltext/S0969-2126(12)00141-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0969212612001414%3Fshowall%3Dtrue

https://www.cell.com/structure/fulltext/S0969-2126(13)00157-3?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0969212613001573%3Fshowall%3Dtrue



Method

The objects to be matched are represented by cubic grids. The center of the grid is at the centroid of the object and grid points are assigned values with regard to their position relative to the object. Grid points within the volume of any virtual atom are considered part of the object (interior or surface); other grid points are assigned “outside the object” values. Distinction between the interior and a surface layer is made for the stationary object A. The resultant grid representations are as follows:

Object A (the stationary object):

![image](https://github.com/user-attachments/assets/66b07116-f245-452d-aede-6286cbf607ec)


Object B (the moving object):

![image](https://github.com/user-attachments/assets/0ab97468-8c9f-47f7-9dbe-417a6313ea45)



The grid representations are correlated using FFT, producing the matrix Cα,β,γ (eq. 1). C holds the correlation scores for shifts of Bl,m,n with respect to Al,m,n by α, β, and γ grid points along three perpendicular axes.

![image](https://github.com/user-attachments/assets/a96ef785-d6f0-4f77-808b-a2db0b9643fc)



By setting Oa to a negative value, s to a positive value and ρ and Ob to zero, the correlation scores are positive when object B is positioned within the volume of object A and overlaps its surface layer. The correlation scores reflect the degree of surface similarity between objects A and B thus higher positive scores indicate more extensive matching of the surfaces and little or no protrusion of object B outside of the volume of object A. The n highest scoring α, β, γ shifts in the current correlation matrix are identified and saved.

To complete a six dimensional search in the rotation-translation space the moving object is rotated to a new orientation (by applying a rotation transformation to the coordinates of the virtual atoms), its grid representation is recalculated for the new orientation, and then a new correlation matrix is calculated and new high scoring positions are identified and saved. The procedure is repeated until the stepwise rotational scan is completed and then all the saved models are sorted by their correlation scores.
