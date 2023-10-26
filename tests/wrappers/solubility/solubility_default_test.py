import json
import os
import requests
import time
import unittest

V2_HOST = os.environ.get("V2_HOST", "0.0.0.0")
V2_PORT = os.environ.get("V2_PORT", "9100")


class SolubilityTest(unittest.TestCase):
    """Test class for Solubility wrapper"""

    @classmethod
    def setUpClass(cls) -> None:
        """This method is run once before all tests in this class."""
        cls.session = requests.Session()
        cls.v2_url = f"http://{V2_HOST}:{V2_PORT}/api/legacy/solubility/batch"

    def get_result(self, task_id: str, timeout: int = 20):
        """Retrieve celery task output"""
        # Try to get result 10 times in 2 sec intervals
        for _ in range(timeout // 2):
            response = self.session.get(f"{self.v2_url}/retrieve?task_id={task_id}")
            result = response.json()
            if result.get("complete"):
                return result
            else:
                if result.get("failed"):
                    self.fail("Celery task failed.")
                else:
                    time.sleep(2)

    def test_1(self):
        case_file = "tests/wrappers/solubility/default_test_case_1.json"
        with open(case_file, "r") as f:
            data = json.load(f)

        # get sync response
        response_sync = self.session.post(
            f"{self.v2_url}/call-sync", json=data
        ).json()

        # get async response
        task_id = self.session.post(
            f"{self.v2_url}/call-async", json=data
        ).json()
        print(task_id)
        time.sleep(10)
        response_async = self.session.get(
            f"{self.v2_url}/retrieve?task_id={task_id}"
        ).json()

        print(response_sync)
        print(response_async)
        # for response in [response_sync, response_async]:
        #     self.assertEqual(response["status_code"], 200)
        #     self.assertIsInstance(response["result"], list)
        #     self.assertIsInstance(response["result"][0], list)
        #     self.assertEqual(len(response["result"][0]), 123)

