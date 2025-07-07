library(data.table)
library(ggplot2)
library(simDAG)

# Fonction pour générer un DAG avec un nombre spécifique de nœuds gaussiens
generate_gaussian_dag <- function(num_vars) {
  # Créer un DAG vide
  dag <- empty_dag()
  
  # Ajouter le nœud racine Bernoulli avec une probabilité de succès de 0.5
  dag <- dag + node("A", type = "rbernoulli", p = 0.5)
  
  # Ajouter les nœuds gaussiens
  for (i in 1:num_vars) {
    node_name <- paste0("Y", i)
    dag <- dag + node(node_name, type = "gaussian", parents = c("A"), betas = c(1), intercept = 0, error = 2)
  }
  
  return(dag)
}

# Fonction pour générer un DAG avec un nombre spécifique de nœuds non gaussiens
generate_non_gaussian_dag <- function(num_vars) {
  # Créer un DAG vide
  dag <- empty_dag()
  
  # Ajouter le nœud racine Bernoulli avec une probabilité de succès de 0.5
  dag <- dag + node("A", type = "rbernoulli", p = 0.5)
  
  # Ajouter les nœuds non gaussiens
  for (i in 1:num_vars) {
    node_name <- paste0("Y", i)
    dag <- dag + node(node_name, type = "gamma", parents = c("A"), betas = c(1), intercept = 0, error = 2)
  }
  return(dag)
}

# Nombre de DAG à générer
num_dags <- 3

# Nombre de variables pour chaque type de DAG (5 avec des variables gaussiennes et 5 avec des variables non gaussiennes)
#num_vars <- round(seq(2, 15, length.out = 5))
num_vars <- round(seq(5, 15, length.out = 3))

# Générer les 5 DAGs avec des variables gaussiennes
dag_gaussian_1 <- generate_gaussian_dag(num_vars[1])
dag_gaussian_2 <- generate_gaussian_dag(num_vars[2])
dag_gaussian_3 <- generate_gaussian_dag(num_vars[3])

# Générer les 5 DAGs avec des variables non gaussiennes
#dags_non_gaussian <- lapply(num_vars, generate_non_gaussian_dag)

# Fusionner les deux listes de DAGs
#all_dags <- c(dag_gaussian, dags_non_gaussian)

# Afficher les DAGs
#print(all_dags)

# data
sim_dat_gauss_1 <- sim_from_dag(dag_gaussian_1, n_sim=10000)
#sim_dat_gamma_1 <- sim_from_dag(dags_non_gaussian[[1]], n_sim=10000)
adj_1 <- dag2matrix(dag_gaussian_1, 
                    include_root_nodes=TRUE, 
                    include_td_nodes=FALSE)
sim_dat_gauss_2 <- sim_from_dag(dag_gaussian_2, n_sim=10000)
sim_dat_gauss_3 <- sim_from_dag(dag_gaussian_3, n_sim=10000)

#1st graph
library(fairadapt)
mod <- fairadapt(Y5 ~ ., 
                 train.data = sim_dat_gauss_1,
                 prot.attr = "A", adj.mat = adj.mat,
                 quant.method=linearQuants)
adapt.train <- adaptedData(mod)

mod2 <- fairadapt(X2 ~ ., 
                  train.data = df,
                  prot.attr = "S", adj.mat = adj.mat,
                  quant.method=rangerQuants)
adapt.train2 <- adaptedData(mod2)

