IFS=","

for args in \
  20,generic,generic_worker \
  2,context_recommender,context_recommender_worker \
  1,evaluate_reactions,evaluate_reactions_worker \
  10,forward,forward_worker \
  1,general_selectivity,general_selectivity_worker \
  1,impurity_predictor,impurity_predictor_worker \
  1,pathway_ranker,pathway_ranker_worker \
  10,retro,retro_worker \
  1,site_selectivity,site_selectivity_worker \
  10,tree_analysis,tree_analysis_worker \
  1,tree_optimizer,tree_optimizer_worker \
  10,tree_search_expand_one,tree_search_expand_one_worker \
  2,tree_search_mcts,tree_search_mcts_worker; \
  do set -- $args;
  CONCURRENCY=$1;
  QUEUE_NAME=$2;
  WORKER_NAME=$3;
  celery -A askcos2_celery worker -c "$CONCURRENCY" -Q "$QUEUE_NAME" -n "$WORKER_NAME"@%h --pool=gevent
  done;
