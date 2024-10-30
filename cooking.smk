vegetables = ['tomatoes', 'onions', 'carots', 'garlic']

rule all:
    input:
        expand("{veggies}_chopped",
            veggies=vegetables)

rule chop_vegetables:
    input:
        "{vegetable}"
    output:
        "{vegetable}_chopped"
    where:
        "frying_pan"
    utensils:
        "knife"
        "chopping_board"
    what:
        """
        Use knife to chop vegetables into small chunks
        """

rule fry_tomatoes:
    input:
        "tomatoes_chopped"
        "oil"
    output:
        "tomatoes_fried"
    where:
        "frying_pan"
    utensils:
        "spatula"
    what:
        """
        Add tomatoes to frying pan and fry for 10mins
        """

rule simmer_sauce:
    input:
        "tomatoes_fried",
        "vegetable_broth"
    output:
        "tomato_sauce"
    where:
        "frying_pan"
    utensils:
        "spatula"
    what:
        """
        Simmer sauce for 20mins
        """

rule taste: 
    input:
        "tomatoes",
        "tomatoes_chopped",
        "tomato_sauce"
    output:
        "taste_summary"
    where:
        "chef"
    utensils:
        "spoon"
    what:
        """
        Taste tomatoes
        """
