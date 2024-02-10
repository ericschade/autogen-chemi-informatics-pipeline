import autogen
from autogen.cache import Cache
import agent_functions
from typing_extensions import Annotated
import asyncio

config_list = autogen.config_list_from_json(env_or_file="OAI_CONFIG_LIST")

llm_config = {
    "config_list": config_list,
}

coder = autogen.AssistantAgent(
    name="chembl_data_engineer_assistant",
    system_message="For coding tasks, only use the functions you have been provided. You have the ability to download data about protein targets and compounds from the chembl database.",
    llm_config=llm_config,
    description="An Agent that is proficient at requesting and saving data from the Chembl web database.",
    human_input_mode="NEVER",
)

plotter = autogen.AssistantAgent(
    name="cheminformatics_plotter",
    system_message="You have the ability to plot chemical properties. Suggest installing dependencies for plotting.",
    llm_config=llm_config,
    description="An Agent that is proficient at visualizing cheminformatics data.",
    human_input_mode="NEVER",
)

# create a UserProxyAgent instance named "user_proxy"
user_proxy_for_running_code = autogen.UserProxyAgent(
    name="user_proxy_for_running_code",
    system_message="A proxy for the user for executing code. When a chembl query is successful, give a short summary of the results returned by the function. When a plot is generated, share the name of the file that is generated.",
    is_termination_msg=lambda x: x.get("content", "") and x.get("content", "").rstrip().endswith("TERMINATE"),
    human_input_mode="NEVER",
    max_consecutive_auto_reply=10,
    code_execution_config={"work_dir": "data_pull/", "use_docker": False},
    description="An agent that should only be used to execute code suggested by assistant agents."
)

user_proxy_for_user_input = autogen.UserProxyAgent(
    name="user_proxy_for_user_input",
    system_message="A proxy for the user for conversation.",
    is_termination_msg=lambda x: x.get("content", "") and x.get("content", "").rstrip().endswith("TERMINATE"),
    human_input_mode="ALWAYS",
    code_execution_config=False,
    description="An agent that should be used when an assistant requests more information from the user."
)

# create an AssistantAgent instance named "cheminformatics_workflow_manager"
cheminformatics_workflow_manager = autogen.AssistantAgent(
    name="cheminformatics_workflow_manager",
    system_message="You are the manager of a cheminformatics drug discovery workflow in which lipinsky descriptors will be calculated for a set of chemical compounds. You help by gathering context about the workflow from the user and efficiently guiding them through the workflow from start to finish, but you do not make tool call suggestions. After data has been saved to a csv by other agents, remind the user that the data is available for them to peruse before moving to the next step in the workflow.",
    llm_config=llm_config,
    description="The manager of a cheminformatics workflow execution. This agent guides the user through the process of data collection and analysis.",
    human_input_mode="NEVER",
)

# define functions according to the function description
autogen.agentchat.register_function(
    agent_functions.download_protein_results,
    caller=coder,
    executor=user_proxy_for_running_code,
    description="Download query results about protein targets from the Chembl web database and save it to a file. The format of the results CSV file is target_query_results_{target_query_string}.csv.",
)

autogen.agentchat.register_function(
    agent_functions.generate_activity_data,
    caller=coder,
    executor=user_proxy_for_running_code,
    description="Download query results about compound activity related to a specific protein target from the Chembl web database and save it to a file. The format of the results CSV file is activity_data_{chembl_id}_{standard_type}.csv.",
)

autogen.agentchat.register_function(
    agent_functions.select_target_from_query_results,
    caller=coder,
    executor=user_proxy_for_running_code,
    description="Get data about a specific protein target from a locally saved csv file of protein target query results.",
)

autogen.agentchat.register_function(
    agent_functions.calculate_lipinksi_descriptors,
    caller=coder,
    executor=user_proxy_for_running_code,
    description="Calculate lipinsky descriptors for a protein target from a locally saved csv file of protein specific activity data. Saves the results to a file with the format {target_id}_{activity_type}_lipinski.csv.",
)

# Add a function for robust group chat termination
@user_proxy_for_running_code.register_for_execution()
@cheminformatics_workflow_manager.register_for_llm()
@coder.register_for_llm(description="terminate the group chat")
def terminate_group_chat(message: Annotated[str, "Message to be sent to the group chat."]) -> str:
    return f"[GROUPCHAT_TERMINATE] {message}"

groupchat = autogen.GroupChat(
    agents=[
        user_proxy_for_running_code,
        user_proxy_for_user_input,
        coder,
        cheminformatics_workflow_manager,
        plotter,
    ],
    messages=[],
    max_round=100
    )

llm_config_manager = llm_config.copy()
llm_config_manager.pop("functions", None)

manager = autogen.GroupChatManager(
    groupchat=groupchat,
    llm_config=llm_config_manager,
    is_termination_msg=lambda x: "GROUPCHAT_TERMINATE" in x.get("content", ""),
    human_input_mode="NEVER",
)

async def main():
    message = """
    Help a user through the following cheminformatics workflow:
    1) Select a protein target.
    2) Using the chosen protein target, get chembl data on chemical compounds that have been screened against the target protein.
    3) Calculate the Lipinski descriptors for the compound.
    4) Perform additional analysis on the compound.
    5) Terminate the chat.
    """

    with Cache.disk():
        await user_proxy_for_running_code.a_initiate_chat(
            manager,
            message=message
        )

if __name__ == "__main__":
    asyncio.run(main())


