import json
import traceback
from fastapi import Depends, HTTPException, Response
from pydantic import BaseModel
from typing import Annotated, Literal
from utils.oauth2 import oauth2_scheme
from utils.registry import get_util_registry
from wrappers import register_wrapper
from wrappers.base import BaseWrapper
from wrappers.tree_analysis.tb_count_analogs import _tb_count_analogs
from wrappers.tree_analysis.tb_pathway_ranking import _tb_pathway_ranking
from wrappers.tree_analysis.tb_pmi_calculation import _tb_pmi_calculation
from wrappers.tree_analysis.tb_reaction_classification \
    import _tb_reaction_classification


class TreeAnalysisInput(BaseModel):
    result_id: str
    task: Literal[
        "pathway_ranking",
        "reaction_classification",
        "pmi_calculation",
        "count_analogs",
    ]

    # PMI calculator and count_analogs option
    index: int = -1

    # Pathway ranking options
    cluster_trees: bool = True
    cluster_method: Literal["hdbscan", "kmeans"] = "hdbscan"
    cluster_min_samples: int = 5
    cluster_min_size: int = 5

    # Count analogs options
    min_plausibility: float = 0.1
    atom_map_backend: Literal["indigo", "rxnmapper", "wln"] = "rxnmapper"


class TreeAnalysisOutput(BaseModel):
    placeholder: str


class TreeAnalysisResponse(BaseModel):
    placeholder: str


@register_wrapper(
    name="tree_analysis_controller",
    input_class=TreeAnalysisInput,
    output_class=TreeAnalysisOutput,
    response_class=TreeAnalysisResponse
)
class TreeAnalysisController(BaseWrapper):
    """Tree Analysis Controller"""
    prefixes = ["tree_analysis/controller", "tree_analysis"]
    backend_wrapper_names = {
        "pathway_ranking": "pathway_ranker",
        "reaction_classification": "reaction_classification",
        "pmi_calculation": "pmi_calculator",
        "count_analogs": "count_analogs"
    }

    def __init__(self):
        pass        # TODO: proper inheritance

    def call_sync(
        self,
        input: TreeAnalysisInput,
        token: Annotated[str, Depends(oauth2_scheme)]
    ) -> Response:
        results_controller = get_util_registry().get_util(
            module="tree_search_results_controller"
        )
        resp = results_controller.retrieve(result_id=input.result_id, token=token)
        if resp.status_code == 404:
            return resp

        try:
            result = json.loads(resp.body)
        except Exception as e:
            traceback.print_exc()
            return resp

        if result["result_type"] == "ipp":
            raise HTTPException(
                status_code=405,
                detail="Analysis of IPP results not supported."
            )

        result_id = input.result_id
        task = input.task

        output = {
            "result_id": result_id,
            "task": task,
            "success": True,
            "error": None,
        }

        try:
            tb_result = result["result"]
            new_revision = result.get("revision", 0) + 1
            if tb_result.get("version") != 2:
                output["success"] = False
                output["error"] = \
                    "Unable to perform analysis tasks on legacy format " \
                    "tree builder results."

                return Response(
                    content=json.dumps(output),
                    status_code=500,
                    media_type="application/json"
                )
        except Exception as e:
            output["success"] = False
            output["error"] = "Unable to load requested result."
            print("Analysis task failed for tree builder result:", str(e))

            return Response(
                content=json.dumps(output),
                status_code=500,
                media_type="application/json"
            )

        if task == "pathway_ranking":
            result, info = _tb_pathway_ranking(
                tb_result=tb_result,
                clustering=input.cluster_trees,
                cluster_method=input.cluster_method,
                min_samples=input.cluster_min_samples,
                min_cluster_size=input.cluster_min_size
            )
        elif task == "reaction_classification":
            result, info = _tb_reaction_classification(tb_result=tb_result)
        elif task == "pmi_calculation":
            result, info = _tb_pmi_calculation(tb_result=tb_result, index=input.index)
        elif task == "count_analogs":
            result, info = _tb_count_analogs(
                tb_result=tb_result,
                index=input.index,
                min_plausibility=input.min_plausibility,
                atom_map_backend=input.atom_map_backend
            )
        else:
            output["success"] = False
            output["error"] = "Unrecognized tree analysis task type."

            return Response(
                content=json.dumps(output),
                status_code=405,
                media_type="application/json"
            )

        output.update(info)
        updated_result = {
            "result": result,
            "revision": new_revision
        }
        if info["success"]:
            resp = results_controller.update(
                result_id=input.result_id,
                updated_result=updated_result,
                token=token
            )
            if not resp.status_code == 200:
                return resp

        return Response(
            content=json.dumps(output),
            status_code=200,
            media_type="application/json"
        )

    async def call_async(
        self,
        input: TreeAnalysisInput,
        priority: int = 0,
        token: Annotated[str, Depends(oauth2_scheme)] = None
    ) -> str:
        from askcos2_celery.tasks import tree_analysis_task
        async_result = tree_analysis_task.apply_async(
            args=(self.name, input.dict(), token), priority=priority)
        task_id = async_result.id

        return task_id

    async def retrieve(self, task_id: str) -> Response | None:
        return await super().retrieve(task_id=task_id)
