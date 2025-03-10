import os
import requests

class Icav2Auth:
    api_key = os.getenv("API_KEY")
    jwt = "jwt"
    base_url = "https://ica.illumina.com/ica/rest"
    
    def __init__(self):
        print(f"Initialising Icav2Auth...")
        
    def get_jwt(self):
        tokens_endpoint = "/api/tokens"
        url = self.base_url + tokens_endpoint
        headers = {
            'accept: application/vnd.illumina.v3+json',
            'X-API-Key: ' + self.api_key
        }
        print(f"Fetching JWT...")
        response = requests.post(
            url,
            data="",
            headers=headers,
        )
        print(f"{response}")
        print(f"{response["token"]}")
        return response
        

if __name__ == '__main__':
    icav2_auth = Icav2Auth()
    icav2_auth.get_jwt()
        