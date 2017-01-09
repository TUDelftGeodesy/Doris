import xml.etree.ElementTree as ET
import os


def init_datastack_xml(datastack_folder):
    xml_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'stack_info_template.xml')
    tree = ET.parse(xml_file)
    settings = tree.getroot()

    settings.find('.project_path').text = datastack_folder
    settings.find('.initialize_flag').text = 'True'

    input = False
    while input == False:
        user_input = raw_input("Enter the path to the orbit files: ")
        if os.path.exists(user_input):
            settings.find('.orbit_path').text = user_input
            input = True
        else:
            print('The path is incorrect, use another path')

    input = False
    while input == False:
        user_input = raw_input("Enter the path your raw datafiles: ")
        if os.path.exists(user_input):
            settings.find('.data_path').text = user_input
            input = True
        else:
            print('The path is incorrect, use another path')

    input = False
    while input == False:
        user_input = raw_input("Do you want to use parallel computing (True/False): ")
        if user_input == 'True' or user_input == 'False':
            settings.find('.parallel').text = user_input
            input = True
        else:
            print('You should use either True or False')

    if user_input == 'True':
        nodes = raw_input("How many cores do you want to use: ")
        settings.find('.nr_of_jobs').text = nodes

    tree.write(open(os.path.join(datastack_folder, 'stack_info.xml'), 'w'))

    return nodes
