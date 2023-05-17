

# sankey plot -------------------------------------------------------------

package <- c("tidyr", "reshape2", "ggforce")
install.packages(package)
library(reshape2)
library(ggforce)
library(tidyr)
dat = melt(pp)

data <- reshape2::melt(Titanic)
head(data)
#  Class    Sex   Age Survived value
# 1   1st   Male Child       No     0
# 2   2nd   Male Child       No     0
# 3   3rd   Male Child       No    35
# 4  Crew   Male Child       No     0
# 5   1st Female Child       No     0
# 6   2nd Female Child       No     0
dat = gather_set_data(dat,1:2)
data <- gather_set_data(data, 1:4)
head(data)
#   Class    Sex   Age Survived value id     x    y
# 1   1st   Male Child       No     0  1 Class  1st
# 2   2nd   Male Child       No     0  2 Class  2nd
# 3   3rd   Male Child       No    35  3 Class  3rd
# 4  Crew   Male Child       No     0  4 Class Crew
# 5   1st Female Child       No     0  5 Class  1st
# 6   2nd Female Child       No     0  6 Class  2nd
ggplot(data, aes(x, id = id, split = y, value = value)) +
  geom_parallel_sets(aes(fill = Sex), alpha = 0.3, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.3, colour = "black", fill = "white" ) + 
  # colour和fill分别修改边框色和填充色
  geom_parallel_sets_labels(colour = 'black', angle = 0)+
  # angle修改方框填充文字的角度
  theme_classic()
