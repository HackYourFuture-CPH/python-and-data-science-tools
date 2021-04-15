import random
# Homework exercises for Week #1
# Write a python script to print your name and age
from typing import Any, Union

my_info = {'name': 'Malek', 'age': '47'}
print(f'My name is: {my_info["name"]} and my age is :{my_info["age"]} year.')

# Create a list of your 5 favorite movies and store it in the variable
my_movie_list = {'name': ['Reader', 'Three Colors', 'Poaching', 'Winter, Spring, Summer, Fall and Winter']}
print(f'{my_movie_list["name"]}')

# Write a Python program to display the first and last colors from the following list.
color_list = ["Red", "Green", "White", "Black"]
print(f'The first color is {color_list[0]} and the last color is {color_list[3]}')

# Write a Python script to add a key to a dictionary


numbers = {0: 10, 1: 20}
numbers.update({'2': 30})
print(numbers)

# Write a Python program to calculate body mass index.
weight = 71
height = 1.73
BMI = weight / height
print(round(BMI, 2))

# Additional Exercises:

# copied!
target_num, guess_num = random.randint(1, 10), 0
while target_num != guess_num:
    guess_num = int(input('Guess a number between 1 and 10 until you get it right : '))
print('Well guessed!')

# Create a list of 5 city names and convert it into tuples.
tuples_test = ['string', 47, 4.7, True]
data_type = list(map(type, tuples_test))
print("The data types are : " + str(data_type))

# Create a list of 5 city names and convert it into tuples.


city_name = ['Copenhagen', 'Tokyo', 'Beijing', 'Moscow', 'Washington dc']
print(tuple(city_name))
# Remove duplicated from the list:


a = [10, 20, 30, 20, 10, 50, 60, 40, 80, 50, 40]
a = list(dict.fromkeys(a))
print(a)


# Accept a string from user and remove the characters which have odd index values of a given string and print them.
def odds(input):
    result = ""
    for i in range(len(input)):
        if i % 2 == 0:
            result: Union[input, Any] = result + input[i]
    return result


print(odds(input("Inter something to get back odds: ")))