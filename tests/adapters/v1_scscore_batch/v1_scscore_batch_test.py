import json
import os
import requests
import time
import unittest

V2_HOST = os.environ.get("V2_HOST", "0.0.0.0")
V2_PORT = os.environ.get("V2_PORT", "9100")


class V1SCScoreBatchTest(unittest.TestCase):
    """Test class for V1 SCScore adapter"""

    @classmethod
    def setUpClass(cls) -> None:
        """This method is run once before all tests in this class."""
        cls.session = requests.Session()
        cls.v1_username = "askcos"
        cls.v1_password = "MML4chem"

        cls.v1_scscore_batch_url = "https://askcos-demo.mit.edu/api/v2/scscore/batch/"
        cls.v1_celery_url = "https://askcos-demo.mit.edu/api/v2/celery/task"
        cls.legacy_adapter_url = f"http://{V2_HOST}:{V2_PORT}/api/legacy/scscore/batch/"
        cls.legacy_celery_url = f"http://{V2_HOST}:{V2_PORT}/api/legacy/celery/task"

    def get_result(self, task_id: str, celery_url: str, mode: str, timeout: int = 20):
        """Retrieve celery task output"""
        # Try to get result 10 times in 2 sec intervals
        for _ in range(timeout // 2):
            if mode == "v1":
                response = self.session.get(
                    url=f"{celery_url}/{task_id}/",
                    auth=(self.v1_username, self.v1_password)
                )
            else:
                response = self.session.get(f"{celery_url}/{task_id}/")

            result = response.json()
            if result.get("complete"):
                return result
            else:
                if result.get("failed"):
                    self.fail("Celery task failed.")
                else:
                    time.sleep(2)

    def test_1(self):
        case_file = "tests/adapters/v1_scscore_batch/v1_scscore_batch_test_case_1.json"
        with open(case_file, "r") as f:
            data = json.load(f)

        results = {}
        for mode, url in [("v1", self.v1_scscore_batch_url),
                          ("legacy", self.legacy_adapter_url)]:
            if mode == "v1":
                response = self.session.post(
                    url=url,
                    json=data,
                    auth=(self.v1_username, self.v1_password)
                )
            else:
                response = self.session.post(url, json=data)

            # Copied from askcos_site/api2/api_test.py
            self.assertEqual(response.status_code, 200)

            # Confirm that request was interpreted correctly
            result = response.json()
            request = result["request"]
            self.assertEqual(request["smiles"], data["smiles"])

            # Test that we got the correct results
            print(result)
            self.assertIsInstance(result["result"], dict)
            results[mode] = result
            print(mode)
            print(results)

        for r1, r2 in zip(results["v1"]["request"]["smiles"], results["legacy"]["request"]["smiles"]):
            self.assertAlmostEqual(results["v1"]["result"][r1], results["legacy"]["result"][r2], places=4)
