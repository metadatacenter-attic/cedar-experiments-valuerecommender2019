# Evaluation utils

import sys
import json
import string
import urllib.parse
from jsonpath_rw import parse
import arm_constants

MISSING_VALUE = 'NA'


def get_test_instances_folder(database, annotated_instances):
    """
    Returns the path to the appropriate test instances folder
    :param database: 
    :param is_training: If False, then it's testing
    :param annotated_instances: 
    :return: 
    """
    if database == arm_constants.BIOSAMPLES_DB.NCBI:
        db = 'NCBI'
    elif database == arm_constants.BIOSAMPLES_DB.EBI:
        db = 'EBI'
    else:
        raise Exception('Invalid database')

    if not annotated_instances:  # free text
        return arm_constants.EVALUATION_TESTING_INSTANCES_BASE_FOLDERS[db]
    else:
        return arm_constants.EVALUATION_TESTING_INSTANCES_ANNOTATED_BASE_FOLDERS[db]


# Returns the instance fields (with their types and values) that are used for the evaluation as a Panda data frame
def get_instance_fields_types_and_values(instance, field_details):
    field_values = {}
    for field in field_details:
        jsonpath_expr = parse(field_details[field]['json_path'])
        matches = jsonpath_expr.find(instance)
        field_values[field] = {}
        field_values[field]['fieldType'] = None
        field_values[field]['fieldValueLabel'] = None
        field_values[field]['fieldValueType'] = None
        if matches[0] is not None and matches[0].value is not None:
            if '@value' in matches[0].value and matches[0].value['@value'] is not None:
                field_values[field]['fieldValueLabel'] = matches[0].value['@value']
            elif 'rdfs:label' in matches[0].value and matches[0].value['rdfs:label'] is not None:
                field_values[field]['fieldValueLabel'] = matches[0].value['rdfs:label']
            if '@id' in matches[0].value and matches[0].value['@id'] is not None:
                field_values[field]['fieldValueType'] = matches[0].value['@id']
            if '@type' in matches[0].value and matches[0].value['@type'] is not None:
                field_type = matches[0].value['@type']
                field_values[field]['fieldType'] = get_original_term_uri(field_type)

    return field_values


# Get the original term uris from BioPortal URIs
def get_original_term_uri(field_type):
    substr = '/classes/'
    if substr in field_type:
        index = field_type.find(substr)
        encoded_term_uri = field_type[index + len(substr):]
        return urllib.parse.unquote(encoded_term_uri)  # return the decoded uri
    else:
        return field_type


def get_populated_fields(field_details, field_values, target_field):
    """
    Returns the populated fields body, given the target field
    :param field_details: 
    :param field_values: 
    :param target_field: 
    :return: 
    """
    populated_fields = []
    for f in field_values:
        if f != target_field:
            if field_values[f] is not None:
                populated_fields.append({'path': field_details[f]['path'], 'value': field_values[f]})
    return populated_fields


# def get_field_path(field_name):
#     return field_paths[field_name]['path']


# Checks if all field values are filled out
def all_filled_out(field_values):
    for f in field_values:
        if field_values[f] is None or len(field_values[f]) == 0:
            return False
    return True


# Extracts the right value from a 'fields_types_and_values' object
def get_field_value(fields_types_and_values, field_name):
    if field_name in fields_types_and_values:
        if 'fieldValueType' in fields_types_and_values[field_name] and fields_types_and_values[field_name]['fieldValueType'] is not None:
            return fields_types_and_values[field_name]['fieldValueType']
        elif 'fieldValueLabel' in fields_types_and_values[field_name] and fields_types_and_values[field_name]['fieldValueLabel'] is not None:
            return fields_types_and_values[field_name]['fieldValueLabel']
        else:
            sys.exit('field value not found')
    else:
        sys.exit('field name not found')


def get_recommended_values(recommendation_results, max_size, annotated):
    if max_size <= 0:
        raise ValueError("max_size must be > 0")
    recommended_values = []
    if 'recommendedValues' in recommendation_results:
        for rv in recommendation_results['recommendedValues']:
            if 'valueType' in rv and rv['valueType'] is not None and annotated:
                recommended_values.append(rv['valueType'])
            elif 'valueLabel' in rv:
                recommended_values.append(rv['valueLabel'])
            else:
                sys.exit('value not found')
              
    else:
        print('Error: recommendedValues not found in recommendation_results')
    result = recommended_values[:max_size]
    return result


def get_recommended_values_as_string(recommended_values):
    if len(recommended_values) == 0:
        return MISSING_VALUE
    elif len(recommended_values) == 1:
        return recommended_values[0]
    else:
        return "|".join(recommended_values)


def is_same_concept(term_uri1, term_uri2, mappings):
    """
    Checks if two term uris have the same meaning. It makes use of a file with mappings
    :param term_uri1: 
    :param term_uri2: 
    """
    if term_uri1 == term_uri2:
        return True
    elif (term_uri2 in mappings and term_uri1 in mappings[term_uri2]) \
            or (term_uri1 in mappings and term_uri2 in mappings[term_uri1]):
        #print('Found two uris for the same concept: ' + term_uri1 + ' = ' + term_uri2)
        return True
    else:
        return False


def get_matching_score(expected_value, value, mappings, normalization=True, extend_with_mappings=False):
    if value == MISSING_VALUE:
        return MISSING_VALUE
    else:
        if extend_with_mappings:
            if is_same_concept(value, expected_value, mappings):
                return 1
            else:
                return 0
        else:
            if normalization:
                value = value.lower()
                expected_value = expected_value.lower()

                # Remove punctuation
                value = value.translate(str.maketrans('', '', string.punctuation))
                expected_value = expected_value.translate(str.maketrans('', '', string.punctuation))

            if expected_value == value:
                return 1
            else:
                return 0


def populated_fields_to_string(populated_fields):
    if len(populated_fields) == 0:
        return 'NA'
    field_value_pairs = []
    for pf in populated_fields:
        field_value_type = ''
        if 'fieldValueType' in pf and pf['fieldValueType'] is not None:
            field_value_type = '[' + pf['fieldValueType'] + ']'
        field_value_pairs.append(pf['fieldPath'] + '=' + field_value_type + pf['fieldValueLabel'])
    result = " | ".join(field_value_pairs)
    return result


def position_of_expected_value(expected_value, values, mappings, extend_with_mappings=False):
    position = 1  # 1-based position
    for value in values:
        if get_matching_score(expected_value, value, mappings, extend_with_mappings=extend_with_mappings) == 1:
            return position
        position += 1
    return 'NA'  # If not found


# Calculates the Reciprocal Rank (https://en.wikipedia.org/wiki/Mean_reciprocal_rank). The MRR will be computed later
def reciprocal_rank(number_of_positions, expected_value, actual_values, mappings, use_na=False,
                    extend_with_mappings=False):
    if use_na and (actual_values == MISSING_VALUE or actual_values is None or len(actual_values) == 0):
        return MISSING_VALUE
    else:
        values = actual_values[0, number_of_positions]
        position = 1
        for value in values:
            if get_matching_score(expected_value, value, mappings, extend_with_mappings) == 1:
                return 1 / float(position)
            position += 1
        return 0


# Note that position is 1-based
def reciprocal_rank_using_position(number_of_positions, position, use_na=False):
    if position is None or position is 'NA':
        if use_na:
            return MISSING_VALUE
        else:
            return 0
    else:
        if position > number_of_positions:
            return 0
        else:
            return 1 / float(position)


def save_to_folder(instance, instance_number, output_path, output_base_file_name):
    """
    Saves an instance to a local folder
    :param instance: 
    :param instance_number: Number used to name the output files
    :param output_path: 
    :param output_base_file_name:
    """
    output_file_path = output_path + "/" + output_base_file_name + "_" + str(instance_number) + '.json'

    with open(output_file_path, 'w') as output_file:
        json.dump(instance, output_file, indent=4)
