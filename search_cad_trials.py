# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "requests",
# ]
# ///

import requests
import json
import sys

def search_trials(condition):
    url = "https://clinicaltrials.gov/api/v2/studies"
    params = {
        "query.cond": condition,
        "filter.overallStatus": "RECRUITING",
        "pageSize": 5
    }
    
    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        data = response.json()
        
        print(f"Total Recruiting {condition} Trials: {data.get('totalCount', 0)}")
        print(f"\n--- TOP RECRUITING TRIALS FOR {condition.upper()} ---")
        for study in data.get('studies', []):
            protocol = study.get('protocolSection', {})
            nct_id = protocol.get('identificationModule', {}).get('nctId')
            title = protocol.get('identificationModule', {}).get('briefTitle')
            phase = protocol.get('designModule', {}).get('phases', ['N/A'])
            
            print(f"[{nct_id}] {title}")
            print(f"Phase: {', '.join(phase)}")
            print("-" * 20)
            
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    condition = sys.argv[1] if len(sys.argv) > 1 else "Coronary Artery Disease"
    search_trials(condition)
