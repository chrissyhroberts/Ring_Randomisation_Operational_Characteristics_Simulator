# Example input data dictionary

`example_ring_sizes.csv` is a synthetic dataset created only to demonstrate the required input structure. It contains no real people, index cases or outbreak records.

| Column | Type | Required | Description |
|---|---|---:|---|
| `r_case_index_id` | character | yes | Unique synthetic ring identifier |
| `direct_contact` | integer | yes | Number of eligible direct contacts; must be greater than zero |
| `contact_of_contact` | integer | yes | Number of eligible contacts-of-contacts; must be zero or greater |

Replace this file with a setting-appropriate empirical distribution before a substantive run.

