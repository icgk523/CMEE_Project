import os
import csv
import re

data_directory = "../Results/TOGA_Output_Apis"
data_sub_directories = [subdir for subdir in os.listdir(data_directory) if
                        os.path.isdir(os.path.join(data_directory, subdir))]

output_dict = {}

for subdir in data_sub_directories:
    if subdir.startswith("TOGA_output_GCA_"):
        toga_name = re.match(r"TOGA_output_(GCA_\d+\.\d+)", subdir).group(1)
        toga_directory = os.path.join(data_directory, subdir)

        # Find the log file in the TOGA directory
        log_files = [file for file in os.listdir(toga_directory) if file.startswith("toga_") and file.endswith(".log")]
        if log_files:
            log_file_path = os.path.join(toga_directory, log_files[0])

            # Read the last 15 lines of the log file
            with open(log_file_path, 'r') as f:
                last_lines = f.readlines()[-15:]

            # Find and extract orthology classification information from the last 15 lines
            orthology_counts = {'one2one': 0, 'one2zero': 0, 'many2one': 0, 'one2many': 0, 'many2many': 0}
            for line in last_lines:
                match = re.match(r"\* (\w+): (\d+)", line.strip())
                if match:
                    classification, count = match.groups()
                    if classification in orthology_counts:
                        orthology_counts[classification] += int(count)

            # Store orthology counts in the output dictionary
            output_dict[toga_name] = orthology_counts

        # Process loss_summ_data.tsv file
        loss_file = os.path.join(toga_directory, "loss_summ_data.tsv")
        if os.path.exists(loss_file):
            with open(loss_file, "r", newline='') as csvfile:
                reader = csv.reader(csvfile, delimiter='\t')
                for row in reader:
                    if row and row[0].startswith("GENE") and len(row) >= 3:
                        value = row[2].strip()
                        output_dict.setdefault(toga_name, {}).setdefault(value, 0)
                        output_dict[toga_name][value] += 1

# Define the path to save the summary table
output_csv_file = "summary_table_hymenoptera.csv"

# Write the summary table to a CSV file
with open(output_csv_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    header = ["Accession", "one2one", "one2zero", "many2one", "one2many", "many2many"]
    extra_columns = sorted(
        set(value for counts in output_dict.values() for value in counts.keys() if value not in header))
    header.extend(extra_columns)
    writer.writerow(header)
    for toga_name, counts in output_dict.items():
        row = [toga_name, counts.get('one2one', 0), counts.get('one2zero', 0), counts.get('many2one', 0),
               counts.get('one2many', 0), counts.get('many2many', 0)]
        row.extend([counts.get(value, 0) for value in extra_columns])
        writer.writerow(row)

print(f"Summary table saved to {output_csv_file}")

