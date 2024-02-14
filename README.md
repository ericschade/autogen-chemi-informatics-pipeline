# AI Tool for chemiinformatics
 
 A repository for analyzing protein/compound activity using a group of Autogen/OpenAi chat agents

![Group Chat Visualization](cheminformatics_group_chat.png)


## Current Issues

Here are some of the current issues we are working on:

- Unwanted conversation termination
    - the groupchat manager is terminating the conversation prematurely, usually just after the plotter agent suggests code but before it is executed.
- GPT4 seems to know some CHEMBL IDs for common protein targets. This means that the workflow manager will sometimes mention the id without asking the user for it. this is problematic because it does not properly allow for user input after each step of the workflow.
- Hallucination of undefined function calls... not sure how to explain it but im seeing many function calls along the lines of "scatter_plot_lipinski_DataVisualization_assistant" or "get_compound_data" when they dont exist
