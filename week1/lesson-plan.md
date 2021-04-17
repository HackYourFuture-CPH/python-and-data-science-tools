# Week 1 : Introduction about course and Python 

## Pre-requisites:
Installation Python and Pycharm

## Lesson Plan:
### 1.	Introduction to Machine learning + application
 [slides](https://github.com/HackYourFuture-CPH/Machine_Learning/blob/main/week1/Python-ML-Week1-Slides.pdf)
 
### 2.	Python installation + Pycharm

Installation Guide : Install [python 3.75](https://www.python.org/downloads/release/python-375/) and [Pycharm community (free) version for practical concepts](https://www.jetbrains.com/pycharm/download/). 

### 3.	Python basic course  

#### Variables
- - - -

We use variables to temporarily store data in the computer’s memory

```
price = 10
rating = 4.9
course_name = ´Python for Begginers´
is_published = True
```

In the example above, <br>

<space>      &emsp; &emsp; <strong>price</strong> is an integer, a whole number without a decimal point;
<br>  &emsp; &emsp; <strong>rating</strong> is a float, a number with a decimal point;
<br>  &emsp; &emsp; <strong>course_name</strong> is a string, a sequence of characters;
<br>  &emsp; &emsp; <strong>is_published</strong> is a boolean. Boolean values can be True or False.

#### Comments
- - - -

We use comments to add notes to our code. Good comments explain the how’s and why’s, not what the code does. That should be reflected in the code itself. Use comments to add reminders to yourself or other developers, or also to explain your assumptions and the reasons you have written code in a certain way.<br>
<br>  &emsp; &emsp; \# This is a comment and it won’t get executed.
<br>  &emsp; &emsp; \# Our comments can be multiple lines.

#### Receiving Input
- - - -

We can receive input from the user by calling the <strong>input()</strong> function.<br>

```
birth_year = int(input(‘Birth year: ‘))
```

The <strong>input()</strong> function always returns data as a string. So, we’re converting the result into an integer by calling the built-in <strong>int()</strong> function.


#### Strings   
- - - -

We can define strings using single (‘ ’) or double (“ ”) quotes.

To define a multi-lined string, we surround our string with triple quotes (“ “ “).

We can get individual characters in a string using square brackets [ ]. <br> 

``` 
course  = ‘Python for Begginers’
course [0]                         \# returns the first character
course [1]                         \# returns the second character
course [-1]                        \# returns the first character from the end
course [-2]                        \# returns the second character from the end
```

We can slice a string using a similar notation:

``` 
course [1:5]
```

The above expression returns all the characters starting from the index position of 1 to 5 (but excluding 5).<br>
The result will be <strong>ytho</strong>

If we leave out the start index, 0 will assumed
If we leave out the end index, the length of the string will be assumed.

We can use formatted strings to dynamically insert values into our strings:

```
name = ‘Mosh’
message = f’Hi, my name is {name } ’

message.upper()                 \# to convert to uppercase
message.lower()                 \# to convert to lowercase
message.title()                 \# to capitalise the first letter of every word
message.find(‘p’)               \# return the index of the first occurrence of p (or -1 not found)
message.replace(‘p’, ‘q’)       \# replaces all the p's for q's
```

To check if a string contains a character or a sequence of characters, we use the <strong> in </strong> operator<br>

```
contains = ‘Python’ in course
```

#### Arithmetic Operations
- - - -

| Operator | Name               | Example  |  Result |
| ----     | -------------- | -----------------------------------| -------- |
| +	| Addition |	x + y |  |
| -	| Subtraction	| x - y |  |
| *	| Multiplication	| x * y | |
| /	| Division	| x / y | returns a float |
| //	| Floor division	| x // y | returns an int |
| %	| Modulo 	| x % y  | returns the remainder of a division|
| %	| Exponentiation | x ** y = x to the power of y <br> or <br> pow(x,y) (pow(base, exponent)) | returns the result of raising the first operand <br>to the power of the second operand.| 


Augmented assignment operator:

```
  x = x + 10<br>
  x += 10
```

Operator precedence
1. Parenthesis
2. Exponentiation
3. Multiplication/Division
4. Addition/Subtraction


#### If Statements     
- - - -

```
if is_hot: 
  print(“hot day”)
elseif is_cold
  print(“cold day”)
else: 
  print(“beautiful day”)
  ```

Logic operators

```
if has_high_income and has_good_credit:
	…
if has_high_income or has_good_credit:
	…
is_day = True
is_night = not is_day
```

#### Comparison operators    
- - - -

| Comparison | Description               | Example  |  Result |
| ----     | ------------------ | ---------| -------- |
| a > b	 |  If the value of a is less than the value of b, then condition becomes true. |	6 > 3 | true  |
| a >= b |  If the value of a is greater than or equal to the value of b, then condition becomes true. |	6 >= 6 | true  |
| a < b  |  If the value of a is less than the value of b, then condition becomes true. |  3 < 6  | true  |
| a <= b |  If the value of a is less than or equal to the b, then condition becomes true. | 6 <= 6 | true  |
| a == b |  If a and b are equal, then the condition becomes true. |6 == 6 | true  |
| a != b |  If a and b are not equal, then condition becomes true. | 6 != 3 | true  |
| a <> b |  If a and b are not equal, then condition becomes true. | 6 <> 3 | true (this operator is deprecated so use != instead) |



#### While loops
- - - -

```
i = 1
while i < 5; 
   print(i)
   i += 1
```

#### For loops     
- - - -

```
for i in range(1, 5):
	print(i)

range(5):         \# generates 0,1,2,3,4
range(1,5):       \# generates 1,2,3,4
range(1,5,2):     \# generates 1,3
```


#### Lists
- - - -
```
numbers = [1, 2, 3, 4, 5]
numbers = [0]                     \#  returns the first item of the list
numbers = [1]                     \#  returns the second item of the list
numbers = [-1]                    \#  returns the first item from the of the list
numbers = [-2]                    \#  returns the second item from the of the list



numbers.append(6)                 \#  adds a single element with the value 6 to the end
numbers.insert(0,6)               \#  adds a single element with the value 6 to the index position 0
numbers.remove(6)                 \#  removes a single and the first element that has the value 6 
numbers.pop()                     \#  removes the last item
numbers.clear()                   \#  removes all the items
numbers.index(8)                  \#  returns the index of the first item with the value 8
numbers.sort()                    \#  returns the list sorted alphabetically 
numbers.reverse()                 \#  returns the list reversed without taking into account the content
numbers.copy()                    \#  returns a copy of the list  
```



#### Tuples     
- - - -

They are like read-only lists. We use them to store a list of items. But once we define a tuple, we cannot add or remove items or change the existing items.

 &emsp; &emsp;  coordinates = (1,2,3)

We can unpack a list or a tuple into separate variables:

 &emsp; &emsp;  x, y, z = coordinates




#### Dictionaries    
- - - -
 
We use dictionaries to store key/value pairs

```
customer =  { 
   “name”: “John Smith”,
   “age”: 30,
   "is_verified”: True
} 
```

We can use strings or numbers to define keys and they should be unique. We can use any type for the values.

```
customer [“name”]                      \# returns “John Smith” 
customer [“type”]                      \# returns an error 
customer.get(“type”, “silver”)         \# returns “silver” 
customer [“name”]  = “new name” 
```


#### Functions    
- - - -

We use functions to break up our code into small chunks. These chunks are easier to read, understand and maintain. If there are bugs, it’s easier to find them in a small chunk than on an entire program. We can also re-use these chunks of code.

```
def great_user(name):
	print(f”Hi {name} ”)

great_user(“John”)
```

<strong>Parameters</strong> are placeholders for the data we can pass to functions. <strong>Arguments</strong> are the actual values we pass.

We have two types of arguments:

* Positional arguments: their position/order matters;
* Keywords arguments: position doesn’t matter and we prefix them with the parameter name.


