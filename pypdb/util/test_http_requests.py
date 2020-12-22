import time
import requests
import unittest
from unittest import mock
import warnings

from pypdb.util import http_requests


class TestHTTPRequests(unittest.TestCase):
    @mock.patch.object(warnings, "warn", autospec=True)
    @mock.patch.object(time, "sleep", autospec=True)
    def test_fails_with_invalid_request(self, mock_sleep, mock_warnings):
        self.assertIsNone(
            http_requests.request_limited(url="http://protein_data_bank.com",
                                          rtype="MAIL"))
        mock_warnings.assert_called_once_with("Request type not recognized")
        self.assertEqual(len(mock_sleep.mock_calls), 0)

    @mock.patch.object(requests, "get", autospec=True)
    @mock.patch.object(time, "sleep", autospec=True)
    def test_get__first_try_success(self, mock_sleep, mock_get):
        mock_response = mock.create_autospec(requests.models.Response)
        mock_response.status_code = 200  # A-OK!
        mock_get.return_value = mock_response

        self.assertEqual(
            http_requests.request_limited(url="http://get_your_proteins.com",
                                          rtype="GET"), mock_response)
        mock_get.assert_called_once_with("http://get_your_proteins.com")
        self.assertEqual(len(mock_sleep.mock_calls), 0)

    @mock.patch.object(requests, "post", autospec=True)
    @mock.patch.object(time, "sleep", autospec=True)
    def test_post__first_try_success(self, mock_sleep, mock_post):
        mock_response = mock.create_autospec(requests.models.Response)
        mock_response.status_code = 200  # A-OK!
        mock_post.return_value = mock_response

        self.assertEqual(
            http_requests.request_limited(url="http://get_your_proteins.com",
                                          rtype="POST"), mock_response)
        mock_post.assert_called_once_with("http://get_your_proteins.com")
        self.assertEqual(len(mock_sleep.mock_calls), 0)

    @mock.patch.object(requests, "get", autospec=True)
    @mock.patch.object(time, "sleep", autospec=True)
    def test_get__succeeds_third_try(self, mock_sleep, mock_get):
        # Busy response
        mock_busy_response = mock.create_autospec(requests.models.Response)
        mock_busy_response.status_code = 429
        # Server Error response
        mock_error_response = mock.create_autospec(requests.models.Response)
        mock_error_response.status_code = 504
        # All good (200)
        mock_ok_response = mock.create_autospec(requests.models.Response)
        mock_ok_response.status_code = 200

        # Mocks `requests.get` to return Busy, then Server Error, then OK
        mock_get.side_effect = [
            mock_busy_response, mock_error_response, mock_ok_response
        ]

        self.assertEqual(
            http_requests.request_limited(url="http://get_your_proteins.com",
                                          rtype="GET"), mock_ok_response)
        self.assertEqual(len(mock_get.mock_calls), 3)
        mock_get.assert_called_with("http://get_your_proteins.com")
        # Should only sleep on being throttled (not server error)
        self.assertEqual(len(mock_sleep.mock_calls), 1)

    @mock.patch.object(warnings, "warn", autospec=True)
    @mock.patch.object(requests, "post", autospec=True)
    @mock.patch.object(time, "sleep", autospec=True)
    def test_post__repeatedly_fails_return_nothing(self, mock_sleep, mock_post,
                                                   mock_warn):
        # Busy response
        mock_busy_response = mock.create_autospec(requests.models.Response)
        mock_busy_response.status_code = 429
        mock_post.return_value = mock_busy_response

        self.assertIsNone(
            http_requests.request_limited(url="http://protein_data_bank.com",
                                          rtype="POST"))
        mock_warn.assert_called_with(
            "Too many failures on requests. Exiting...")

        self.assertEqual(len(mock_post.mock_calls), 4)
        mock_post.assert_called_with("http://protein_data_bank.com")
        self.assertEqual(len(mock_sleep.mock_calls), 4)


if __name__ == '__main__':
    unittest.main()
