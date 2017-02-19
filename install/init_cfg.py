import xml.etree.ElementTree as ET
import os

def init_cfg():
    xml_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'grs_config.xml')
    tree = ET.parse(xml_file)
    settings = tree.getroot()

    settings.find('.source_path').text = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    input = False
    while input == False:
        user_input = raw_input("Enter the path to doris: ")
        if os.path.exists(user_input) and user_input.endswith('doris'):
            settings.find('.doris_path').text = user_input
            input = True
        else:
            print('The path is incorrect, use another path')

    input = False
    while input == False:
        user_input = raw_input("Enter the path to cpxfiddle: ")
        if os.path.exists(user_input) and user_input.endswith('cpxfiddle'):
            settings.find('.cpxfiddle_path').text = user_input
            input = True
        else:
            print('The path is incorrect, use another path')

    tree.write(open(xml_file, 'w'))

# Actually execute the code...
if __name__ == "__main__":

    # Initialize...
    init_cfg()
