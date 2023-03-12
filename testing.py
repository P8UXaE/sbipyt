
# # This is for saving the model generated
# # save the model
# torch.save(gnn_model.state_dict(), 'gnn_model.pt')


# create a new instance of the GNN model
new_gnn_model = GNN(input_dim, hidden_dim, output_dim)

# load the saved state dictionary into the new model
new_gnn_model.load_state_dict(torch.load('gnn_model.pt'))

# use the new model for prediction
new_output = new_gnn_model(test_feature_matrix, test_adj_matrix)
