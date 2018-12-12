from datetime import datetime
import shutil
from os import listdir, remove, walk, mkdir
from os.path import isfile, isdir, join, dirname, exists
import pickle
import glob

"""
EXAMPLE USAGE:
bed_directory_path = "/Users/saureous/data/bed_alleles_copy"
output_file_path = "output/bed_output/concatenated.bed"

file_manager = FileManager()
file_paths = file_manager.get_file_paths_from_directory(directory_path=bed_directory_path)
file_manager.concatenate_files(file_paths=file_paths, output_file_path=output_file_path)
file_manager.delete_files(file_paths=file_paths)
"""


class FileManager:
    """
    Does simple file operations like concatenation, fetching a list of paths for files in a directory, deleting files
    """
    @staticmethod
    def concatenate_files(file_paths, output_file_path):
        """
        Concatenate files (efficiently) given a list of file paths to a single file
        :param file_paths: List of file path
        :param output_file_path: Output file path name
        :return: None
        """
        with open(output_file_path, 'wb') as out_file:
            for file_path in file_paths:
                with open(file_path, 'rb') as in_file:
                    # 100MB per writing chunk to avoid reading big file into memory.
                    shutil.copyfileobj(in_file, out_file, 1024*1024*100)

    @staticmethod
    def get_file_paths_from_directory(directory_path):
        """
        Returns all paths of files given a directory path
        :param directory_path: Path to the directory
        :return: file_paths: A list of paths of files
        """
        file_paths = [join(directory_path, file) for file in listdir(directory_path) if isfile(join(directory_path, file))]
        return file_paths

    @staticmethod
    def get_subdirectories(directory_path):
        """
        Collect the absolute paths of all top-level directories contained in the specified directory
        :param directory_path: Path to a directory containing subdirectories
        :return: dir_paths: list of paths
        """
        dir_paths = [join(directory_path, file) for file in listdir(directory_path) if isdir(join(directory_path, file))]
        return dir_paths

    @staticmethod
    def get_all_file_paths_by_type(parent_directory_path, file_extension, sort=True):
        """
        Given a parent directory, iterate all files within, and return those that end in the extension provided by user.
        File paths returned in sorted order by default.
        :param parent_directory_path:
        :param file_extension:
        :param sort:
        :return:
        """
        all_files = set()

        for root, dirs, files in walk(parent_directory_path):
            for file in files:
                if file.endswith(file_extension):
                    all_files.add(join(root,file))

        all_files = list(all_files)

        if sort:
            all_files.sort()

        return all_files

    @staticmethod
    def ensure_directory_exists(directory_path, i=0):
        """
        Recursively test directories in a directory path and generate missing directories as needed
        :param directory_path:
        :return:
        """
        if i > 3:
            print("WARNING: generating subdirectories of depth %d, please verify path is correct"%i)

        if not exists(directory_path):
            try:
                mkdir(directory_path)

            except FileNotFoundError:
                FileManager.ensure_directory_exists(dirname(directory_path), i=i+1)

                if not exists(directory_path):
                    mkdir(directory_path)

    @staticmethod
    def delete_files(file_paths):
        """
        Deletes files given in file paths
        :param file_paths: List of file paths
        :return: None
        """
        for file_path in file_paths:
            remove(file_path)

    @staticmethod
    def get_datetime_string():
        now = datetime.now()
        now = [now.year, now.month, now.day, now.hour, now.minute, now.second, now.microsecond]
        datetime_string = "_".join(list(map(str, now)))

        return datetime_string

    @staticmethod
    def save_object_pickle(output_dir, filename, object):
        array_file_extension = ".pkl"

        # ensure chromosomal directory exists
        if not exists(output_dir):
            FileManager.ensure_directory_exists(output_dir)

        output_path_prefix = join(output_dir, filename)

        output_path = output_path_prefix + array_file_extension

        with open(output_path, 'wb') as output:
            pickle.dump(object, output, pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    test_directory = "./test/directory/path/a/b/"
    FileManager.ensure_directory_exists(test_directory)
