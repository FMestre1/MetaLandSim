\name{spom}
\alias{spom}
\title{
Stochastic Patch Occupancy Model
}
\description{
This function predicts the occupancy status of each patch in a landscape in the time step t+1, based on the occupancy information on time step t.
}
\usage{
spom(sp, kern, conn, colnz, ext, param_df, 
     beta1 = NULL, b = 1, c1 = NULL, c2 = NULL, 
     z = NULL, R = NULL, succ="none", max_age=1)
}
\arguments{
  \item{sp}{
Landscape with species occupancy, object of class 'metapopulation'. 
}
  \item{kern}{
'op1' or 'op2'. Dispersal kernel. See details.
}
  \item{conn}{
'op1' or 'op2'. Connectivity function. See details.
}
  \item{colnz}{
'op1', 'op2' or 'op3'. Colonization function. See details.
}
  \item{ext}{
'op1', 'op2' or 'op3'. Extinction function. See details.
}
  \item{param_df}{
Parameter data frame delivered by \code{\link{parameter.estimate}}, including:
\itemize{
 \item alpha - Parameter relating extinction with distance.
 \item y - Parameter y in the colonization probability.  
 \item e - Parameter defining the extinction probability in a patch of unit area.
 \item x - Parameter scaling extinction risk with patch area.
 }
}
  \item{beta1}{
Parameter affecting long distance dispersal probability (if the Kern='op2').
}
  \item{b}{
Parameter scaling emigration with patch area (if conn='op1' or 'op2'). By default set to 1.
}
  \item{c1}{
Parameter scaling immigration with the focal patch area (if conn='op2'). 
}
  \item{c2}{
Parameter c in the option 3 of the colonization probability (if colnz='op3').
}
  \item{z}{
Parameter giving the strength of the Allee effect (if colnz='op3').
}
  \item{R}{
Parameter giving the strength of the Rescue effect (if ext='op3').
}
  \item{succ}{
Set the preference of the species for patch successional stage: 'none', 'early', 'mid' and 'late'.
}
  \item{max_age}{
Default value set to 1. This argument should not be changed by the user. It is used only when the function runs inside others.
}
}
\details{
In order to visualize which parameter combination is valid for each option, please refer to the following table (alpha, x, y and e are delivered by \code{\link{parameter.estimate}}, as a data frame):

\tabular{ccccccccccc}{
\strong{parameter}  \tab kern_1 \tab kern_2 \tab conn_1 \tab conn_2 \tab colnz_1 \tab colnz_2 \tab colnz_3 \tab ext_1 \tab ext_2\tab ext_3\cr
\strong{alpha} \tab x \tab x \tab \tab \tab \tab \tab \tab \tab \tab \cr
\strong{x} \tab \tab \tab \tab \tab \tab \tab \tab x \tab x \tab x \cr
\strong{y} \tab \tab \tab \tab \tab x \tab x \tab \tab \tab \tab \cr
\strong{e} \tab \tab \tab \tab \tab \tab \tab \tab x \tab x \tab x \cr
\strong{beta1} \tab \tab x \tab \tab \tab \tab \tab \tab \tab \tab \cr
\strong{b} \tab \tab \tab x \tab x \tab \tab \tab \tab \tab \tab \cr
\strong{c1} \tab \tab \tab \tab x \tab \tab \tab \tab \tab \tab \cr
\strong{c2} \tab \tab \tab \tab \tab \tab \tab x \tab \tab \tab \cr
\strong{z} \tab \tab \tab \tab \tab \tab \tab x \tab \tab \tab \cr
\strong{R} \tab \tab \tab \tab \tab \tab \tab \tab \tab \tab x \cr
}

A Stochastic Patch Occupancy Model (SPOM) is a type of model which models the occupancy status of the species on habitat patches as a Markov chain (Moilanen, 2004). These models are a good compromise between capturing sufficient biological detail and  being easy to parametrize with occupancy data. With SPOMs it is possible to predict the probability of extinction or colonization of every patch in a landscape, given the current occupancy state of all the patches (Etienne et al. 2004).

\emph{Dispersal Kernel}


Option 1 (Hanski, 1994 and 1999)
\deqn{D(D_{ij},\alpha) = exp(-\alpha.d_{ij})}{D(Dij,\alpha) = exp(-\alpha.dij)}

Option 2 (Shaw, 1995)
\deqn{D(D_{ij},\alpha,\beta) = \frac{1}{1+\alpha.d_{ij}^\beta}}{D(Dij,\alpha)=1/(1+\alpha.dij^\beta)}

where \emph{dij is the distance between patches i and j.}

\itemize{
 \item Option 1 - Negative exponential. Earlier studies (until the end of the 1990) frequently used this type of thin-tailed kernels (Nathan et al. 2012).
 \item Option 2 - Fat-tailed kernel. The shape of the dispersal kernel is highly significant only when the metapopulation consists of several moderately small patch clusters, which are relatively far from each other. In this kind of a system, a patch cluster may go extinct, and long-distance dispersal will be important in determining the recolonization probability of the empty cluster (Shaw, 1995 and Moilanen, 2004). This type of fat-tailed kernels has become more frequent in recent works (Nathan et al. 2012). For \deqn{\beta=2} this is the Cauchy distribution. 
}

\emph{Connectivity}


Option 1 (Moilanen, 2004)
\deqn{S_i=\sum pj.D(d_{ij},\alpha).A_j^b}{Si=sum(pj.D(dij,alpha).Aj^b)}

Option 2 (Moilanen and Nieminen, 2002)
\deqn{S_i=A_i^c \sum p_j.D(d_{ij},\alpha).A_j^b}{Si=Ai^c sum(pj.D(dij,alpha).Aj^b)}

where \emph{Ai and Aj are the areas of patches i(focal patch) and j(other patches), respectively; dij is the distance between patches i and j and pj is the occupation status (0/1) of patch j}\cr

\itemize{
 \item Option 1 - In the version of Hanski (1994), de kernel is the negative exponential (option 1) and b is set to 1. In this more flexible version, the parameter b scales emigration with patch area (Moilanen, 2004).
 \item Option 2 - In Moilanen & Nieminen (2002) the kernel is the negative exponential (option 1). This metric considers the value of the focal patch's area, which was found to provide better results by Moilanen & Nieminen (2002), being less sensitive to errors in the estimation of a. Parameters b and c scale, respectively emigration and immigration, as a function of patch area (focal patch in the case of c). See 'note'.
}

\emph{Colonization function}

Option 1 (Hanski, 1994, 1999)
\deqn{C_i=\frac{S_i^2}{S_i^2+y^2}}{Ci=Si^2/(Si^2+y^2)}

Option 2 (Moilanen, 2004)
\deqn{C_i=1-exp(-y.S_i)}{Ci=1-exp(-y.Si)}

Option 3 (Ovaskainen, 2002)
\deqn{C_i=\frac{S_i^z}{S_i^z+\frac{1}{c}}}{Ci=Si^z/(Si^z+1/c)}


where \emph{Si is connectivity.}


\itemize{
 \item Option 1 - It's the first version of the colonization probability, it includes Allee effect (however the strength of this effect cannot be modified) Hanski (1994). Colonization probability is defined as a sigmoid function of the connectivity of patch i. 
 \item Option 2 - This option assumes that immigrating individuals originate colonization events independently, therefore, with no Allee effect. Adequate for species (plants) with passive dispersal (Moilanen, 2004). 
 \item Option 3 - Here, as in option 1, the colonization probability is defined as a sigmoid function of the connectivity of patch i, and the user can change the strength of the Allee effect, by changing the parameter z, with values >1 reflecting the presence of this effect (Ovaskainen, 2002). In the original version of the IFM (option 1) Hanski (1994) assumed a relatively strong Allee effect (z=2). Parameter c describes the species ability to colonize (Ovaskainen & Hanski, 2001 and Ovaskainen ,2002).
}


\emph{Extinction function}

Option 1 (Hanski, 1994, 1999)
\deqn{E_i=min(1,\frac{e}{A_i^x})}{Ei=min(1,e/Ai^x)}

Option 2 (Hanski and Ovaskainen, 2000 and Ovaskainen and Hanski, 2002)
\deqn{E_i=1-exp(\frac{-e}{A_i^x})}{Ei=1-exp(-e/Ai^x)}

Option 3 (Ovaskainen, 2002)
\deqn{E_i=min[1,\frac{e}{A_i^x}.(1-C_i)^R]}{Ei=min(1,((e/Ai^x).(1-Ci)^R))}

where \emph{Ai is the area of the focal patch and Ci is the colonization probability of the focal patch.}

\itemize{
 \item Option 1 - Original version developed by Hanski (1994).
 \item Option 2 - Used e.g. in the spatially realistic Levins model (Hanski & Ovaskainen, 2000 and Ovaskainen & Hanski, 2002). Parameter x scales extinction probability with patch area.
 \item Option 3 - Same as option 1, but considering the Rescue effect (with the strength of this effect being given by R). If R=0 there is no Rescue effect, however, if R>0, the Rescue effect grows exponentially with the probability of not being colonized. In the original version of this function Hanski (1994) assumed R=1.
}
Here, parameter \emph{x} defines de degree to which the extinction rate is sensitive to the patch area. If \emph{x>1}, with the increase of \emph{Ai} the extinction rate rapidly approximates zero. The populations in the larger patches becomes almost impossible to extinguish. However, if \emph{x} is small the extinction rate decreases slower with increasing \emph{Ai}. 
}
\value{
Delivers a list similar to the class 'metapopulation' but with two additional columns in the data frame nodes.characteristics: 'species2'(which is the occupation in the next time step) and turn (turnover between occupancies).
}
\references{
Etienne, R. S., ter Braak, C. J., and Vos, C. C. (2004). Application of stochastic patch occupancy models to real metapopulations. In Hanski, I. and Gaggiotti, O.E. (Eds.) Ecology, Genetics, and Evolution of Metapopulations. Elsevier Academic Press. 696 pp.

Hanski, I. (1994). A practical model of metapopulation dynamics. Journal of Animal Ecology, 63: 151-162.

Hanski, I. (1999). \emph{Metapopulation Ecology}. Oxford University Press. 313 pp.

Hanski, I., Alho, J., and Moilanen, A. (2000). Estimating the parameters of survival and migration of individuals in metapopulations. Ecology, 81(1), 239-251.

Hanski, I., and Ovaskainen, O. (2000). The metapopulation capacity of a fragmented landscape. Nature, 404: 755-758.

Moilanen, A. (2004). SPOMSIM: software for stochastic patch occupancy models of metapopulation dynamics. Ecological Modelling, 179(4), 533-550.

Moilanen, A., and Nieminen, M. (2002). Simple connectivity measures in spatial ecology. Ecology, 83(4): 1131-1145.

Nathan, R., Klein, E., Robledo-Arnuncio, J.J. and Revilla, E. (2012). Dispersal kernels: review. in Clobert, J., Baguette, M., Benton, T. and Bullock, J.M. (Eds.) Dispersal Ecology and Evolution. Oxford University Press. Oxford, UK. 462 pp.  

Ovaskainen, O. (2002). The effective size of a metapopulation living in a heterogeneous patch network. The American Naturalist: 160(5), 612-628.

Ovaskainen, O. and Hanski, I. (2001). Spatially structured metapopulation models: global and local assessment of metapopulation capacity. Theoretical Population Biology, 60(4), 281-302.

Ovaskainen, O., and Hanski, I. (2002). Transient dynamics in metapopulation response to perturbation. Theoretical Population Biology, 61(3): 285-295.

Ovaskainen, O. and Hanski, I. (2004). Metapopulation dynamics in highly fragmented landscapes. In Hanski, I. & Gaggiotti, O.E. (Eds.) Ecology, Genetics, and Evolution of Metapopulations. Elsevier Academic Press. 696 pp.

Shaw, M.W., (1995). Simulation of population expansion and spatial pattern when individual dispersal distributions do not decline exponentially with distance. Proc. R. Soc. London B: 259, 243-248.
}
\author{
Frederico Mestre and Fernando Canovas
}
\note{
Future versions of the package should include the virtual migration model (Hanski et al. 2000), which allows the estimation of migration related parameters (relevant to the option 2 of connectivity).
}
\seealso{
\code{\link{species.graph}}, \code{\link{simulate_graph}}, \code{\link{iterate.graph}}
}
\examples{

data(occ.landscape)
data(param1)

#Simulating the occupation in the next time step:

landscape2 <- spom(sp=occ.landscape,
			kern="op1",
			conn="op1",
			colnz="op1",
			ext="op1",
			param_df=param1,
			beta1=NULL,
			b=1,
			c1=NULL,
			c2=NULL,
			z=NULL,
			R=NULL,
			succ="none"
			)

#The output has two new columns in the data frame nodes.characteristics: species2 
#(occupation in the next time step) and turn (turnover - change of occupation status, 
#1 if changed and 0 if not).:

head(landscape2)

#         x         y      areas    radius cluster    colour nneighbour
#1 718.5011 228.47190 0.05741039 13.518245       1 #FF0000FF   91.80452
#2 494.3624  73.29165 0.08755563 16.694257       1 #FF0000FF   98.98432
#3 809.2326 245.90046 0.09384384 17.283351       1 #FF0000FF  166.68205
#4 638.8057 149.35122 0.08858989 16.792569       1 #FF0000FF   82.60306
#5 874.2010  19.78104 0.03621793 10.737097       1 #FF0000FF   92.26625
#6 605.3937  70.34944 0.03066018  9.878987       1 #FF0000FF  131.22261
#  ID species species2 turn
#1  1       1        1    0
#2  2       0        1    1
#3  3       1        1    0
#4  4       0        0    0
#5  5       0        1    1
#6  6       1        1    0

}
