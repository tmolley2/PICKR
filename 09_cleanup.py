#09_cleanup.py
import shutil
import os

def delete_folder(folder_path):
    try:
        if os.path.exists(folder_path):
            shutil.rmtree(folder_path)  # Deletes the folder and its contents
            print(f"Folder '{folder_path}' has been deleted successfully.")
        else:
            print(f"Folder '{folder_path}' does not exist.")
    except Exception as e:
        print(f"An error occurred while deleting the folder: {e}")
        
def move_folder_contents(source_folder, destination_folder):
    # Ensure source and destination folders exist
    if not os.path.exists(source_folder):
        print(f"Source folder '{source_folder}' does not exist.")
        return
    if not os.path.exists(destination_folder):
        os.makedirs(destination_folder)  # Create destination folder if it doesn't exist

    # Iterate through items in the source folder
    for item in os.listdir(source_folder):
        source_path = os.path.join(source_folder, item)
        destination_path = os.path.join(destination_folder, item)

        try:
            # Move files or directories
            shutil.move(source_path, destination_path)
            print(f"Moved: {source_path} -> {destination_path}")
        except Exception as e:
            print(f"Error moving {source_path}: {e}")

#move newly run csvs into done data
source = "data/new"
destination = "data"
move_folder_contents(source, destination)

# delete temp files
folder_to_delete = "temp"
delete_folder(folder_to_delete)

